"""A small, dependency-free local web interface for HCRProbeForge.

The browser form is intentionally a front end, not a second implementation:
single-target submissions call ``core.main`` and multi-target submissions call
the same Python batch runner used by the cross-platform helper commands.
"""

from __future__ import annotations

import argparse
import errno
import html
import ipaddress
import json
import mimetypes
import os
import platform
import re
import secrets
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import threading
import webbrowser
from http.cookies import CookieError, SimpleCookie
from email.parser import BytesParser
from email.policy import default
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Callable
from urllib.parse import parse_qs, quote, unquote, urlparse, urlsplit
from uuid import uuid4

from . import batch, core, references


JOBS: dict[str, dict[str, object]] = {}
RUN_LOCK = threading.Lock()
JOBS_LOCK = threading.Lock()
LIFECYCLE_LOCK = threading.Lock()
ACTIVE_JOB_STATUSES = frozenset({"queued", "running", "cancelling"})


WEB_TOKEN_HEADER = "X-HCRProbeForge-Token"
WEB_TOKEN_COOKIE = "HCRProbeForge-Token"
DEFAULT_WEB_PORT = 8766
WEB_PORT_ATTEMPTS = 10


ProgressCallback = Callable[[dict[str, object]], None]


LOGO_SVG = r'''<svg class="logo-mark" viewBox="0 0 500 120" preserveAspectRatio="xMinYMid meet" role="img" aria-labelledby="hcrpf-logo-title hcrpf-logo-desc" xmlns="http://www.w3.org/2000/svg">
<title id="hcrpf-logo-title">HCRProbeForge</title>
<desc id="hcrpf-logo-desc">HCRProbeForge logo showing two probe halves hybridizing to a transcript</desc>
<defs>
  <linearGradient id="hcrpf-logo-bg" x1="0" y1="0" x2="1" y2="1">
    <stop offset="0" stop-color="#f8f7fc"/>
    <stop offset="1" stop-color="#efecf8"/>
  </linearGradient>
</defs>
<circle cx="96" cy="60" r="46" fill="url(#hcrpf-logo-bg)" stroke="#7366a9" stroke-width="4.5"/>
<line x1="68" y1="73" x2="124" y2="73" stroke="#34324b" stroke-width="4.5" stroke-linecap="round"/>
<line x1="76" y1="65" x2="76" y2="70" stroke="#d56b85" stroke-opacity=".5" stroke-width="2" stroke-linecap="round"/>
<line x1="82" y1="65" x2="82" y2="70" stroke="#d56b85" stroke-opacity=".5" stroke-width="2" stroke-linecap="round"/>
<line x1="88" y1="65" x2="88" y2="70" stroke="#d56b85" stroke-opacity=".5" stroke-width="2" stroke-linecap="round"/>
<line x1="104" y1="65" x2="104" y2="70" stroke="#7366a9" stroke-opacity=".5" stroke-width="2" stroke-linecap="round"/>
<line x1="110" y1="65" x2="110" y2="70" stroke="#7366a9" stroke-opacity=".5" stroke-width="2" stroke-linecap="round"/>
<line x1="116" y1="65" x2="116" y2="70" stroke="#7366a9" stroke-opacity=".5" stroke-width="2" stroke-linecap="round"/>
<path d="M72 63h17c6 0 6-16 0-16h-13" fill="none" stroke="#d56b85" stroke-width="4.5" stroke-linecap="round" stroke-linejoin="round"/>
<path d="M120 63h-17c-6 0-6-16 0-16h13" fill="none" stroke="#7366a9" stroke-width="4.5" stroke-linecap="round" stroke-linejoin="round"/>
<text x="152.96" y="79.25" fill="#34324b" font-family="Helvetica-Bold, Helvetica Neue, Helvetica, Arial, sans-serif" font-size="38" font-weight="700">HCR</text>
<text x="236.28" y="79.25" fill="#7366a9" font-family="Helvetica-Bold, Helvetica Neue, Helvetica, Arial, sans-serif" font-size="38" font-weight="700">Probe</text>
<text x="344.96" y="79.25" fill="#d56b85" font-family="Helvetica-Bold, Helvetica Neue, Helvetica, Arial, sans-serif" font-size="38" font-weight="700">Forge</text>
</svg>'''


FAVICON_SVG = r'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="44 8 104 104" role="img" aria-labelledby="hcrpf-favicon-title">
<title id="hcrpf-favicon-title">HCRProbeForge</title>
<defs>
  <linearGradient id="hcrpf-favicon-bg" x1="0" y1="0" x2="1" y2="1">
    <stop offset="0" stop-color="#f8f7fc"/>
    <stop offset="1" stop-color="#efecf8"/>
  </linearGradient>
</defs>
<circle cx="96" cy="60" r="46" fill="url(#hcrpf-favicon-bg)" stroke="#7366a9" stroke-width="4.5"/>
<line x1="68" y1="73" x2="124" y2="73" stroke="#34324b" stroke-width="4.5" stroke-linecap="round"/>
<line x1="76" y1="65" x2="76" y2="70" stroke="#d56b85" stroke-opacity=".5" stroke-width="2" stroke-linecap="round"/>
<line x1="82" y1="65" x2="82" y2="70" stroke="#d56b85" stroke-opacity=".5" stroke-width="2" stroke-linecap="round"/>
<line x1="88" y1="65" x2="88" y2="70" stroke="#d56b85" stroke-opacity=".5" stroke-width="2" stroke-linecap="round"/>
<line x1="104" y1="65" x2="104" y2="70" stroke="#7366a9" stroke-opacity=".5" stroke-width="2" stroke-linecap="round"/>
<line x1="110" y1="65" x2="110" y2="70" stroke="#7366a9" stroke-opacity=".5" stroke-width="2" stroke-linecap="round"/>
<line x1="116" y1="65" x2="116" y2="70" stroke="#7366a9" stroke-opacity=".5" stroke-width="2" stroke-linecap="round"/>
<path d="M72 63h17c6 0 6-16 0-16h-13" fill="none" stroke="#d56b85" stroke-width="4.5" stroke-linecap="round" stroke-linejoin="round"/>
<path d="M120 63h-17c-6 0-6-16 0-16h13" fill="none" stroke="#7366a9" stroke-width="4.5" stroke-linecap="round" stroke-linejoin="round"/>
</svg>'''

FAVICON_TAG = '<link rel="icon" type="image/svg+xml" href="/favicon.svg">'


STYLE = r"""
.check[hidden]{display:none!important}
/* Pastel interface accents shared with the companion HCR applications. */
.workflow div:nth-child(1){background:#f0ebff!important}.workflow div:nth-child(2){background:#e8f7f2!important}.workflow div:nth-child(3){background:#eef4fb!important}.workflow div:nth-child(4){background:#fff0e8!important}.muted{color:var(--muted)}
:root{--ink:#3c3850;--muted:#716a82;--accent:#7968ad;--accent2:#8876bd;--bg:#f8f6fc;--card:#fffdfd;--line:#ebe5f1;--mint:#dff6f0;--sand:#fff1df;--rose:#f9dce8;--lav:#e8e0ff;--peach:#ffe5da}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif}main{max-width:1120px;margin:34px auto;padding:0 22px}.hero{display:flex;justify-content:space-between;gap:24px;align-items:center;margin-bottom:20px}.eyebrow{text-transform:uppercase;letter-spacing:.12em;color:var(--accent);font-size:12px;font-weight:800}.hero h1{font-size:40px;line-height:1.05;margin:8px 0}.subtitle{max-width:760px;color:var(--muted);font-size:17px;line-height:1.45}.version{background:var(--mint);color:var(--accent);font-weight:800;border-radius:999px;padding:7px 11px;font-size:13px;white-space:nowrap}.workflow{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:20px 0}.workflow div{padding:16px;border:1px solid var(--line);border-radius:15px;background:#fff}.workflow b{display:block;margin-bottom:5px}.workflow small{color:var(--muted);line-height:1.35}.card{background:var(--card);border:1px solid var(--line);border-radius:18px;padding:26px;box-shadow:0 10px 32px #5c4c7a12!important;margin:18px 0}.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px;align-items:start}.full{grid-column:1/-1}label{display:block;font-weight:700;margin-bottom:7px}.field-label{display:flex;align-items:center;gap:7px}input,select,textarea{width:100%;padding:11px 12px;border:1px solid #c7d6d1;border-radius:10px;background:#fff;color:var(--ink);font:inherit}textarea{min-height:88px;resize:vertical}input:focus,select:focus,textarea:focus{outline:3px solid #d8eeea;border-color:var(--accent)}select[multiple]{min-height:110px}.help{display:inline-grid;place-items:center;width:18px;height:18px;border-radius:50%;background:#e8e0ff!important;color:#62558b!important;font-size:12px;cursor:help;position:relative}.help:hover::after,.help:focus::after{content:attr(data-help);position:absolute;z-index:4;width:min(360px,75vw);margin:28px 0 0 -20px;padding:12px 14px;background:#3c3850;color:#fff;border-radius:9px;font-size:13px;line-height:1.45;box-shadow:0 8px 20px #5c4c7a35}.note{background:var(--sand);color:#66563c;padding:13px 15px;border-radius:11px;line-height:1.45}.tip{font-size:13px;color:var(--muted);margin-top:5px;line-height:1.4}.mode-note{display:none;background:var(--mint);padding:12px 14px;border-radius:10px;color:#356258}.secondary{background:#fff;color:var(--accent);border:1px solid #a9c7c0;box-shadow:none}.actions{display:flex;gap:12px;align-items:center;margin-top:20px}button{background:linear-gradient(135deg,var(--accent),var(--accent2));border:0;border-radius:11px;color:#fff;padding:13px 20px;font:800 16px inherit;cursor:pointer;box-shadow:0 7px 18px #236b6830}button:hover{filter:brightness(1.05)}summary{cursor:pointer;color:var(--accent);font-weight:800}.inspection{display:none;background:var(--mint);border-radius:10px;padding:11px 13px;color:#356258;white-space:pre-wrap}.error{background:var(--rose);border-left:5px solid #ca6c83;border-radius:10px;padding:16px;white-space:pre-wrap}.console{background:#17242b;color:#d9eee8;padding:16px;border-radius:11px;white-space:pre-wrap;overflow:auto;max-height:460px;font:13px/1.45 ui-monospace,SFMono-Regular,Menlo,monospace}.file-list{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px;margin:0;padding:0;list-style:none}.file-list li{border:1px solid var(--line);border-radius:11px;padding:11px 13px;background:#fbfdfc;min-width:0}.file-list a{color:var(--accent);font-weight:700;text-decoration:none;display:block;overflow-wrap:anywhere;word-break:break-word;line-height:1.3}.file-list a:hover{text-decoration:underline}.file-kind{display:block;color:var(--muted);font-size:12px;margin-top:5px;overflow-wrap:anywhere;word-break:break-word;line-height:1.3}.path{display:block;background:#eee8f5!important;padding:10px 12px;border-radius:9px;overflow:auto;white-space:nowrap;font:13px ui-monospace,SFMono-Regular,Menlo,monospace}.progress-track{height:16px;background:#e9e3f7!important;border-radius:999px;overflow:hidden}.progress-bar{height:100%;width:0;background:linear-gradient(90deg,var(--accent),var(--accent2));transition:width .25s ease}.progress-bar.indeterminate{width:45%;animation:progress-slide 1.35s infinite ease-in-out}@keyframes progress-slide{0%{transform:translateX(-110%)}100%{transform:translateX(245%)}}.status-line{display:flex;justify-content:space-between;gap:12px;margin:10px 0;color:var(--muted)}.map-gallery{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:18px}.map-card{margin:0;border:1px solid var(--line);border-radius:13px;padding:12px;background:#fff!important}.map-card img{display:block;width:100%;height:auto;border-radius:8px}.map-card figcaption{font-size:13px;color:var(--muted);margin-top:8px;word-break:break-word}.result-actions{display:flex;flex-wrap:wrap;gap:10px;align-items:center}.button-link{display:inline-block;background:linear-gradient(135deg,var(--accent),var(--accent2));border-radius:11px;color:#fff;padding:11px 15px;font-weight:800;text-decoration:none}.button-link.secondary{background:#fff;color:var(--accent)}.small-note{font-size:13px;color:var(--muted);line-height:1.45}.cancel-button{background:#fff0f4;color:#bd5b78;border:1px solid #e4a3b5;box-shadow:none}.cancel-button:hover{background:#f9dce8}@media(max-width:760px){.grid,.workflow{grid-template-columns:1fr}.full{grid-column:auto}.hero{display:block}.logo-mark{width:min(360px,100%);height:auto}.version{display:inline-block;margin-top:8px}.file-list{grid-template-columns:1fr}.map-gallery{grid-template-columns:1fr}}
.workflow-tabs-card{background:#fbfaff;border:1px solid var(--line);border-radius:14px;padding:15px 16px}


.workflow-tab.active{background:var(--accent);color:#fff;border-color:var(--accent)}
.mode-select{position:absolute;width:1px;height:1px;overflow:hidden;clip:rect(0 0 0 0);white-space:nowrap}
.guide-link{display:inline;padding:0;border:0;background:none;color:var(--accent);font:700 17px inherit;box-shadow:none;text-decoration:underline;vertical-align:baseline}
.guide-link:hover{filter:none;color:#bd5b78}
.welcome-backdrop{position:fixed;inset:0;z-index:20;display:grid;place-items:center;padding:22px;background:#3c385066;backdrop-filter:blur(3px)}
.welcome-backdrop.hidden{display:none}
.species-modal-backdrop{position:fixed;inset:0;z-index:25;display:grid;place-items:center;padding:22px;background:#3c385066;backdrop-filter:blur(3px)}
.species-modal-backdrop.hidden{display:none}
.species-modal-card{position:relative;width:min(720px,100%);max-height:min(860px,calc(100vh - 44px));overflow:auto;background:var(--card);border:1px solid var(--line);border-radius:20px;padding:30px 32px;box-shadow:0 20px 70px #3c385055}.species-modal-card .help:hover::after,.species-modal-card .help:focus::after{content:none}.species-delete-card{position:relative;width:min(560px,100%);max-height:min(720px,calc(100vh - 44px));overflow:auto;background:var(--card);border:1px solid var(--line);border-radius:20px;padding:30px 32px;box-shadow:0 20px 70px #3c385055}.species-delete-card h2{margin:0 34px 8px 0;font-size:28px}.species-delete-card p{color:var(--muted);line-height:1.5}.delete-confirm{display:flex;align-items:flex-start;gap:12px;padding:14px 0;font-weight:700;line-height:1.45}.delete-confirm input{width:18px;height:18px;flex:0 0 18px;margin:2px 0 0}.delete-message{display:none;margin:4px 0 14px}
.floating-help{position:fixed;z-index:100;width:min(360px,calc(100vw - 24px));padding:12px 14px;background:#3c3850;color:#fff;border-radius:9px;font-size:13px;line-height:1.45;box-shadow:0 8px 20px #5c4c7a35;pointer-events:none}
.species-modal-card h2{margin:0 34px 8px 0;font-size:30px}.species-modal-card p{color:var(--muted);line-height:1.5}.species-modal-card .required{color:#bd5b78}.species-modal-card .optional,.species-modal-card .conditional{font-weight:500;color:var(--muted)}.species-modal-message{display:none;margin-top:14px}.delete-assets{padding:12px 14px;background:var(--lav);border-radius:10px;overflow-wrap:anywhere;word-break:break-word}
.welcome-card{position:relative;width:min(640px,100%);max-height:min(780px,calc(100vh - 44px));overflow:auto;background:var(--card);border:1px solid var(--line);border-radius:20px;padding:30px 32px;box-shadow:0 20px 70px #3c385055;font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif}
.welcome-card h2{margin:0 34px 8px 0;font-size:30px}.welcome-card p{color:var(--muted);line-height:1.5}.welcome-card ol{margin:18px 0 22px;padding-left:24px;color:var(--ink)}.welcome-card li{padding:7px 0;line-height:1.45}.welcome-card li::marker{color:var(--accent);font-weight:800}.welcome-card code{background:var(--lav);padding:2px 5px;border-radius:5px}.welcome-close{position:absolute;top:12px;right:14px;background:none;color:var(--muted);box-shadow:none;padding:4px 9px;font-size:28px;line-height:1}.welcome-close:hover{background:var(--lav);filter:none}.welcome-actions{display:flex;justify-content:flex-end;gap:10px}.welcome-actions button{min-width:120px}.index-status{display:inline-block;margin-left:8px;padding:4px 9px;border-radius:999px;font-size:12px;font-weight:800;background:var(--sand);color:#66563c}.index-status.ready{background:var(--mint);color:#356258}
.brand{align-items:center;display:grid;grid-template-columns:minmax(390px,520px) minmax(0,1fr);gap:18px;min-width:0}.logo-mark{width:min(520px,100%);height:auto;display:block;align-self:center}.header-copy{min-width:0;padding-top:0}.header-copy .subtitle{max-width:650px}.grid>div{min-width:0}.grid>div:not(.full)>label:not(.check){min-height:28px}.grid>#channel_input{grid-column:auto}.grid>#auto_input{grid-column:1/-1;padding-top:0}.check{display:flex!important;align-items:center!important;justify-content:flex-start;gap:12px!important;min-height:46px;padding:0 4px!important;line-height:1.35;min-width:0}.check input{width:18px!important;height:18px!important;margin:0!important;flex:0 0 18px!important;vertical-align:middle}.check>label{display:inline-flex!important;align-items:center;min-width:0;margin:0!important;line-height:1.35;white-space:normal}.check>.help{flex:0 0 18px;align-self:center;margin:0}.full.row{display:flex!important;flex-wrap:wrap;justify-content:center;align-items:center;gap:12px;grid-column:1/-1}.full.row>#inspect{flex:0 1 720px;width:min(100%,720px);margin:0 auto}.full.row>.inspection{flex:1 1 100%;width:100%;margin:0}.row{display:grid;grid-template-columns:minmax(0,1fr) minmax(250px,1fr);gap:9px;align-items:stretch}.row>input{min-width:0}.row>input[type=file]{height:46px;padding:7px 10px}.advanced-grid,.grid#advanced_grid{align-items:start}#advanced_grid[data-mode="design"]>div:not([data-workflows~="design"]),#advanced_grid[data-mode="manifest"]>div:not([data-workflows~="manifest"]),#advanced_grid[data-mode="plot"]>div:not([data-workflows~="plot"]),#advanced_grid[data-mode="qc"]>div:not([data-workflows~="qc"]),#advanced_grid[data-mode="index"]>div,#advanced_grid>div[hidden]{display:none!important}#advanced_grid .check{align-self:start}@media(max-width:900px){.brand{grid-template-columns:minmax(300px,420px) minmax(0,1fr);gap:18px}.logo-mark{width:min(420px,100%)}}@media(max-width:760px){.brand{display:block}.logo-mark{width:min(520px,100%);margin-bottom:14px}.row{grid-template-columns:1fr}.row>input[type=file]{height:auto}.grid>div:not(.full)>label:not(.check){min-height:0}}
.workflow-tabs{display:grid;grid-template-columns:repeat(5,minmax(0,1fr));align-items:stretch;gap:8px;margin-top:9px}.workflow-tab{background:#fff;color:var(--accent);border:1px solid #c9bfe5;box-shadow:none;padding:10px 14px;font-size:14px;width:100%;min-width:0;min-height:44px}.species-picker{display:flex;align-items:center;gap:8px}.species-picker select{flex:1;min-width:0}.species-remove{display:grid;place-items:center;width:42px;height:46px;margin:0;padding:0;font-size:25px;line-height:1;border-radius:10px}.advanced-grid>.check{display:flex!important;align-items:center!important;align-self:start!important;gap:12px!important;min-height:46px;padding:0 4px!important;margin:0!important}.advanced-grid>.check>input{flex:0 0 18px}.advanced-grid>.check>label{flex:0 1 auto}.advanced-grid>.check>.help{flex:0 0 18px;align-self:center}.advanced-grid>.check .help{align-self:center}#index_input .grid>.check{display:grid!important;grid-template-columns:18px minmax(0,1fr) minmax(230px,1fr);align-items:center!important;align-self:start;gap:12px!important;min-height:84px;padding:38px 4px 0!important}.index-folder-button{width:100%;min-width:0;min-height:52px;padding:13px 24px;font-size:17px;line-height:1.2;white-space:nowrap;margin:0}.index-folder-button:hover{filter:none;background:var(--lav)}.index-folder-button:focus{outline:3px solid #d8eeea}@media(max-width:1000px){.workflow-tabs{grid-template-columns:repeat(3,minmax(0,1fr))}}@media(max-width:760px){.workflow-tabs{grid-template-columns:1fr}.advanced-grid>.check{min-height:46px;padding:0 4px!important;transform:none!important}#index_input .grid>.check{display:flex!important;min-height:46px;padding:10px 4px 0!important;flex-wrap:wrap}.index-folder-button{width:100%;margin:0}}
.target-type-panel{grid-column:1/-1;background:#f6f3ff;border:1px solid #e5ddf7;border-radius:13px;padding:15px 16px}.target-type-panel .target-type-row{display:grid;grid-template-columns:minmax(0,1fr) minmax(260px,340px);gap:16px;align-items:end}.target-type-panel .target-type-row label{margin:0}.target-type-panel .target-type-details{margin-top:12px;padding-top:12px;border-top:1px solid #e5ddf7}.target-type-panel .target-type-details[hidden]{display:none}.target-type-panel .target-type-details label{margin-bottom:6px}.target-type-panel .target-type-details input{max-width:380px}.target-type-panel .tip{margin-top:5px}.target-type-panel .target-status{margin-left:7px}.target-type-panel .target-type-copy{max-width:760px}.target-type-panel select{background:#fff}.target-status.ready{background:var(--mint);color:#356258}.target-status.pending{background:var(--sand);color:#66563c}@media(max-width:760px){.target-type-panel .target-type-row{grid-template-columns:1fr}.target-type-panel select{margin-top:2px}}
"""


FORM = r"""
<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">__FAVICON__<title>HCRProbeForge</title><style>__STYLE__</style></head>
<body><main><section class="hero"><div class="brand">__LOGO__<div class="header-copy"><div class="eyebrow">HCR v3 split-initiator probe design</div><div class="subtitle">Choose an organism and workflow, provide your sequence or gene list, and HCRProbeForge will design, check, and save the results on this computer. <button type="button" class="guide-link" id="open_welcome">How to start</button></div></div></div><span class="version">v__VERSION__</span></section>
<section class="workflow"><div><b>Resolve</b><small>Find the requested transcript from NCBI or validate the FASTA you provide.</small></div><div><b>Design</b><small>Generate channel-specific candidates and screen them against the selected genome index.</small></div><div><b>Curate</b><small>Apply sequence filters, oligo QC, fallback tiers, and non-overlapping probe selection.</small></div><div><b>Export</b><small>Save selected pairs, order files, maps, QC workbooks, reports, and reproducibility metadata.</small></div></section>
<form class="card" method="post" action="/run" enctype="multipart/form-data"><div class="grid">
<div class="full workflow-tabs-card"><label class="field-label">Workflow <span class="help" tabindex="0" aria-label="Workflow help" data-help="Choose what you want HCRProbeForge to do. The page will show only the inputs and settings needed for that workflow.">?</span></label><div class="workflow-tabs" role="tablist" aria-label="HCRProbeForge workflows"><button type="button" class="workflow-tab active" data-mode="design">Design one target</button><button type="button" class="workflow-tab" data-mode="manifest">Design a gene list</button><button type="button" class="workflow-tab" data-mode="plot">Plot a probe table</button><button type="button" class="workflow-tab" data-mode="qc">QC an oligo table</button><button type="button" class="workflow-tab" data-mode="index">Build a genome index</button></div><select id="mode" name="mode" class="mode-select"><option value="design">Design one target</option><option value="manifest">Design a gene/accession list</option><option value="plot">Plot an existing probe table</option><option value="qc">QC an oligo table</option><option value="index">Build a genome index</option></select></div>
<div id="reference_input" class="full"><div class="grid"><div id="species_input"><label class="field-label">Species <span class="help" tabindex="0" aria-label="Species help" data-help="Select the organism you are studying. HCRProbeForge uses it to find transcripts and to choose the matching Bowtie2 genome reference.">?</span><span id="index_status" class="index-status" aria-live="polite">Checking index…</span></label><select id="species" name="species"><option data-organism="Xenopus tropicalis" data-assembly-name="UCB_Xtro_10.0" data-assembly-accession="GCF_000004195.4" data-index-ready="__XTR_INDEX_READY__" data-annotation-ready="__XTR_ANNOTATION_READY__" value="xtr">Xenopus tropicalis</option><option data-organism="Xenopus laevis" data-assembly-name="Xenopus_laevis_v10.1" data-assembly-accession="GCF_017654675.1" data-index-ready="__XLA_INDEX_READY__" data-annotation-ready="__XLA_ANNOTATION_READY__" value="xla">Xenopus laevis</option><option data-organism="Danio rerio" data-assembly-name="GRCz12ab" data-assembly-accession="GCF_052040795.1" data-index-ready="__ZEBRAFISH_INDEX_READY__" data-annotation-ready="__ZEBRAFISH_ANNOTATION_READY__" value="zebrafish">Zebrafish (Danio rerio)</option><option data-organism="Mus musculus" data-assembly-name="GRCm39" data-assembly-accession="GCF_000001635.27" data-index-ready="__MOUSE_INDEX_READY__" data-annotation-ready="__MOUSE_ANNOTATION_READY__" value="mouse">Mouse (Mus musculus)</option><option data-organism="Gallus gallus" data-assembly-name="bGalGal1.mat.broiler.GRCg7b" data-assembly-accession="GCF_016699485.2" data-index-ready="__CHICKEN_INDEX_READY__" data-annotation-ready="__CHICKEN_ANNOTATION_READY__" value="chicken">Chicken (Gallus gallus)</option><option data-organism="Homo sapiens" data-assembly-name="GRCh38.p14" data-assembly-accession="GCF_000001405.40" data-index-ready="__HUMAN_INDEX_READY__" data-annotation-ready="__HUMAN_ANNOTATION_READY__" value="human">Human (Homo sapiens)</option>__CUSTOM_SPECIES_OPTIONS__</select><div class="tip">Choose a preset to use its registered index or build its current reference. Use the index workflow for another assembly.</div></div><div id="organism_input"><label class="field-label">NCBI organism <span class="help" tabindex="0" aria-label="NCBI organism help" data-help="This is the organism name sent to NCBI. It is filled automatically from the species choice so transcript searches stay matched to the selected species.">?</span></label><input id="organism" name="organism" value="Xenopus tropicalis" readonly></div></div></div>
<div id="channel_input"><label class="field-label">HCR channel <span class="help" tabindex="0" aria-label="Channel help" data-help="Choose one channel for a smaller run, or All channels to design the target separately with B1, B2, B3, B4, and B5 initiators.">?</span></label><select name="channel"><option value="ALL">All channels (B1-B5)</option><option value="B1">B1</option><option value="B2">B2</option><option value="B3">B3</option><option value="B4">B4</option><option value="B5">B5</option></select></div>
<div id="design_inputs"><label class="field-label">Gene symbol <span class="help" tabindex="0" aria-label="Gene symbol help" data-help="Type the gene symbol you want to design. Leave the accession blank if you want HCRProbeForge to choose a transcript automatically.">?</span></label><input name="gene" placeholder="e.g. sox9 or gnrhr2/nmi"><div class="tip">With only a gene symbol, HCRProbeForge searches NCBI and applies the transcript policy in Advanced settings.</div></div>
<div id="accession_input"><label class="field-label">RefSeq transcript accession <span class="help" tabindex="0" aria-label="Accession help" data-help="Enter a RefSeq RNA accession when you want to design one exact transcript. A version such as NM_001016853.2 is best for reproducible work.">?</span></label><input name="accession" placeholder="e.g. NM_001016853.2"><div class="tip">A versioned accession identifies the exact sequence used for the design.</div></div>
<div id="fasta_input" class="full"><label class="field-label">FASTA file or local FASTA path <span class="help" tabindex="0" aria-label="FASTA help" data-help="Use this when you already have the transcript sequence. Provide one FASTA record, and do not fill in Gene symbol or accession at the same time.">?</span></label><div class="row"><input name="fasta_path" placeholder="Optional local path"><input name="fasta_upload" type="file" accept=".fa,.fasta,.fna,.txt"></div></div>
<div id="target_type_input" class="target-type-panel"><div class="target-type-row"><div class="target-type-copy"><label class="field-label">Design target <span id="premrna_annotation_status" class="index-status target-status pending" hidden aria-live="polite"></span><span class="help" tabindex="0" aria-label="Design target help" data-help="Mature transcript uses processed RNA. Pre-mRNA · intronic regions uses the full genomic transcript but keeps each probe pair inside selected introns. Pre-mRNA · whole genomic transcript allows probes in exons, introns, and across their boundaries.">?</span></label><div id="target_type_tip" class="tip">Mature transcript: probes are designed on the processed, spliced RNA. Choose a Pre-mRNA mode when you need unspliced nuclear RNA; those modes require a matching genomic annotation database.</div></div><select id="target_type" name="target_type"><option value="mature">Mature transcript</option><option value="pre-mrna">Pre-mRNA · intronic regions</option><option value="pre-mrna-whole">Pre-mRNA · whole genomic transcript</option></select></div><div id="premrna_details" class="target-type-details" hidden><label>Introns to target <span class="optional">(optional)</span> <span class="help" tabindex="0" aria-label="Pre-mRNA intron help" data-help="Leave blank to use all annotated introns. Enter 1-based intron numbers separated by commas, such as 1,3,5, to restrict probes to those introns.">?</span></label><input name="premrna_introns" placeholder="All introns, or e.g. 1,3,5"><div class="tip">Leave blank to include every annotated intron. To restrict the design, enter 1,2 or 1,3,5. Probe pairs must fit completely inside the selected introns; exons and exon–intron boundary-spanning pairs are excluded.</div></div></div>
<div id="manifest_input" class="full" style="display:none"><label class="field-label">Gene manifest or list <span class="help" tabindex="0" aria-label="Manifest help" data-help="Upload a spreadsheet or text file with one target per row. Use gene_symbol for the gene, accession when you need one exact transcript, and an optional channel column (B1-B5) when different rows belong to different HCR channels. Rows without a channel use the HCR channel setting above.">?</span></label><div class="row"><input name="manifest_path" placeholder="Optional local path"><input name="manifest_upload" type="file" accept=".tsv,.csv,.xlsx,.xlsm,.txt"></div><div class="tip">Accepted formats are CSV, TSV, XLSX, XLSM, and text. Headers such as <code>gene_symbol</code>, <code>accession</code>, and optional <code>channel</code> are recognized. Species mismatches are reported as skips; unresolved genes and accessions are recorded as failures.</div></div>
<div id="existing_input" class="full" style="display:none"><label class="field-label">Existing input file <span class="help" tabindex="0" aria-label="Existing file help" data-help="Use a probe table for plotting, or an IDT/selected-pairs table for QC. The selected workflow determines how the file is read and what is written.">?</span></label><div class="row"><input name="existing_path" placeholder="Optional local path"><input name="existing_upload" type="file" accept=".tsv,.csv,.xlsx,.xlsm,.txt"></div></div>
<div id="index_input" class="full" style="display:none"><div class="note"><strong>Build a genome index:</strong> select a species above and the preset will supply the compatible NCBI assembly. <span id="index_preset_reference">The selected preset supplies its compatible assembly.</span> For an NCBI build, HCRProbeForge downloads and validates the genome FASTA and matching genomic annotation, builds the Bowtie2 index, creates a compact annotation lookup, and registers the reference. The uncompressed reference files are retained for fast Pre-mRNA designs; duplicate <code>.gz</code> files are removed. A custom/local build can provide both a genome FASTA and its matching GFF3 annotation.</div><div class="grid" style="margin-top:16px"><div><label>NCBI assembly accession <span class="optional">(optional)</span> <span class="help" tabindex="0" aria-label="Assembly accession help" data-help="For a preset with an NCBI assembly, leave this blank to use the preset. For a custom species, provide a versioned GCF_ or GCA_ accession so both genome and annotation are downloaded and checked before indexing.">?</span></label><input name="index_assembly_accession" placeholder="Optional, versioned GCF_ or GCA_ accession"><div class="tip">Leave blank to use the selected species preset.</div></div><div><label>Expected assembly name <span class="optional">(optional)</span> <span class="help" tabindex="0" aria-label="Assembly name help" data-help="If supplied with an NCBI accession, the returned Assembly name must match exactly. This prevents building the wrong reference.">?</span></label><input name="index_assembly_name" placeholder="Optional, e.g. GRCz12ab"><div class="tip">If supplied, the NCBI assembly name must match this value.</div></div><div><label>Index alias <span class="optional">(optional)</span> <span class="help" tabindex="0" aria-label="Index alias help" data-help="The stable alias used by HCRProbeDesign and later design commands. Leave blank to use the species alias for the preset build.">?</span></label><input name="index_alias" placeholder="Optional, e.g. octopus_v1"><div class="tip">Use a new alias to keep an alternate assembly beside the preset index.</div></div><div id="index_local_fasta_input" class="full"><label>Local genome FASTA <span class="optional">(optional)</span> <span class="help" tabindex="0" aria-label="Local genome FASTA help" data-help="Upload a genome FASTA when the assembly is not available from NCBI or when you want to build from your own reference. A saved custom-species FASTA is used automatically when this field is blank.">?</span></label><div class="row"><input name="index_fasta_path" placeholder="Optional local path"><input name="index_fasta_upload" type="file" accept=".fa,.fasta,.fna,.gz"></div></div><div id="index_local_annotation_input" class="full"><label>Matching genomic annotation GFF3 <span class="optional">(recommended for local Pre-mRNA use)</span> <span class="help" tabindex="0" aria-label="Genomic annotation help" data-help="For a local FASTA build, provide the GFF3 from the same assembly. NCBI builds download and prepare this automatically.">?</span></label><div class="row"><input name="index_annotation_path" placeholder="Optional local GFF3 path"><input name="index_annotation_upload" type="file" accept=".gff,.gff3,.gz"></div></div><div><label>Build threads</label><input name="index_threads" type="number" min="1" value="4"><div class="tip">More threads can shorten the build if your computer has them available.</div></div><div class="check"><input name="index_force" type="checkbox"><label style="margin:0">Rebuild an existing index</label></div></div></div>
<div id="output_input"><label class="field-label">Results/project folder <span class="help" tabindex="0" aria-label="Results folder help" data-help="Choose the project folder where HCRProbeForge should create its runs and reusable cache subfolders. Run, cache, plot, QC, and manifest outputs are nested below the folder you enter; the webapp keeps a small index-metadata copy beside the project root.">?</span></label><input name="outdir" value="hcr_results"><div class="tip">Run outputs are stored below <code>runs/</code>; reusable transcript data is stored below <code>cache/</code>.</div></div>
<div id="folder_input"><label class="field-label">Optional workflow folder name <span class="help" tabindex="0" aria-label="Workflow folder help" data-help="Give this run a short name you will recognize later, for example phox2b_minimal. If blank, the input filename or manifest name is used.">?</span></label><input name="folder_name" placeholder="e.g. phox2b_minimal"><div class="tip">The name is placed inside the selected species and workflow folder.</div></div>
<div class="full row"><button class="secondary" type="button" id="inspect">Inspect input</button><div class="inspection" id="inspection"></div></div>
<div id="auto_input" class="check" data-workflows="design manifest"><input id="auto" name="auto_curate" type="checkbox" checked><label for="auto" style="margin:0">Smart auto-curation</label><span class="help" tabindex="0" aria-label="Auto-curation help" data-help="Leave this on for the usual design. If the first pass does not produce enough usable pairs, HCRProbeForge tries its predefined fallback settings and records which pass succeeded.">?</span></div>
<div id="design_note" class="full"><p class="note"><strong>Before designing:</strong> Normal runs need a genome index for the selected species. If the badge says the index isn't registered, run <strong>Build a genome index</strong> first. Turning off genome masking skips this requirement but isn't recommended, since probes won't be screened for off-target matches.</p></div>
<div class="full"><details id="advanced_settings"><summary>Advanced settings</summary><p class="note" style="margin-top:16px"><strong>How these settings work:</strong> the basic fields choose the input and workflow. These optional controls are the same controls available in the command line. They change how transcripts are selected, how candidate probes are filtered and curated, how plots are drawn, or how QC is reported. If you are unsure, keep the defaults; the settings shown here change with the selected workflow.</p><div id="advanced_grid" class="grid" data-mode="design" style="margin-top:16px">
<div data-workflows="design manifest"><label class="field-label">Transcript selection policy <span class="help" tabindex="0" aria-label="Transcript policy help" data-help="For a gene-only request, auto chooses one well-supported linked RefSeq transcript. Choose longest to prioritize length, or require-accession when you want to enter the transcript yourself.">?</span></label><select name="transcript_policy"><option value="auto">auto</option><option value="longest">longest</option><option value="interactive">interactive (CLI only)</option><option value="require-accession">require-accession</option></select></div>
<div class="check" data-workflows="design manifest"><input name="all_transcripts" type="checkbox"><label style="margin:0">Design all linked transcripts</label><span class="help" tabindex="0" aria-label="All transcripts help" data-help="Use this only with a gene symbol and leave Accession blank. HCRProbeForge designs every linked RefSeq RNA transcript instead of choosing one representative transcript.">?</span></div>
<div data-workflows="design manifest"><label class="field-label">Explicit Bowtie2 index prefix <span class="help" tabindex="0" aria-label="Index prefix help" data-help="Normally the selected species supplies the index. Enter a full Bowtie2 prefix here only when you want to use a different registered or local index.">?</span></label><input name="index" placeholder="Optional path; overrides genome alias"></div>
<div data-workflows="design manifest"><label>NCBI email</label><input name="email" placeholder="NCBI_EMAIL or email address"></div>
<div data-workflows="design manifest"><label>NCBI API key</label><input name="api_key" type="password" placeholder="NCBI_API_KEY or key"></div>
<div data-workflows="design manifest"><label class="field-label">Tile size <span class="help" tabindex="0" aria-label="Tile size help" data-help="Length of the target sequence used for each probe pair. The default is 52 nucleotides; shorter tiles are available to smart fallback tiers.">?</span></label><input name="tile_size" type="number" min="1" value="52"></div>
<div data-workflows="design manifest"><label class="field-label">Minimum GC (%) <span class="help" tabindex="0" aria-label="Minimum GC help" data-help="Lowest GC percentage accepted for the target tile in this design pass.">?</span></label><input name="min_gc" type="number" step="any" value="45"></div>
<div data-workflows="design manifest"><label class="field-label">Maximum GC (%) <span class="help" tabindex="0" aria-label="Maximum GC help" data-help="Highest GC percentage accepted for the target tile in this design pass.">?</span></label><input name="max_gc" type="number" step="any" value="55"></div>
<div data-workflows="design manifest"><label>Minimum Gibbs</label><input name="min_gibbs" type="number" step="any" value="-70"></div>
<div data-workflows="design manifest"><label>Maximum Gibbs</label><input name="max_gibbs" type="number" step="any" value="-50"></div>
<div data-workflows="design manifest"><label>Target Gibbs</label><input name="target_gibbs" type="number" step="any" value="-60"></div>
<div data-workflows="design manifest"><label class="field-label">Max C/G-run mismatches <span class="help" tabindex="0" aria-label="C/G run help" data-help="Controls the C/G-rich run filter. A larger value is more permissive; it is not a Bowtie alignment setting.">?</span></label><input name="max_run_mismatches" type="number" min="0" value="2"></div>
<div data-workflows="design manifest"><label class="field-label">Initial max probes <span class="help" tabindex="0" aria-label="Initial max probes help" data-help="How many candidates the initial design pass asks the engine to return. Leave blank to use the automatic value.">?</span></label><input name="max_probes" type="number" min="1" placeholder="Dynamic default"></div>
<div data-workflows="design manifest"><label class="field-label">Maximum genomic hits <span class="help" tabindex="0" aria-label="Genomic hits help" data-help="A candidate is retained only when its genomic alignment count is within this limit. The default of 1 is the strictest specificity setting.">?</span></label><input name="num_hits_allowed" type="number" min="1" value="1"></div>
<div class="check" data-workflows="design manifest"><input name="dtm_filter" type="checkbox"><label style="margin:0">Enable dTm filter</label><span class="help" tabindex="0" aria-label="dTm filter help" data-help="Apply the maximum dTm value to the two probe arms. Leave it off to preserve the design engine's normal behavior.">?</span></div>
<div data-workflows="design manifest"><label class="field-label">Maximum dTm <span class="help" tabindex="0" aria-label="Maximum dTm help" data-help="Maximum allowed difference between the two arm melting temperatures when the dTm filter is enabled.">?</span></label><input name="dtm_max" type="number" min="0" step="any" value="5"></div>
<div class="check" data-workflows="design manifest"><input name="no_genomemask" type="checkbox"><label style="margin:0">Disable genome masking</label><span class="help" tabindex="0" aria-label="Genome masking help" data-help="Skip Bowtie2 specificity screening. Use this only when you deliberately accept that genomic uniqueness was not checked.">?</span></div>
<div data-workflows="design manifest plot"><label class="field-label">Plot theme <span class="help" tabindex="0" aria-label="Plot theme help" data-help="Pastel uses the publication-style colored map; minimal uses a simpler low-decoration layout. This changes appearance only.">?</span></label><select name="plot_theme"><option value="pastel">pastel</option><option value="minimal">minimal</option></select></div>
<div data-workflows="design manifest plot"><label class="field-label">Plot color by <span class="help" tabindex="0" aria-label="Plot color help" data-help="Choose what colors the probe blocks: gc for GC percentage, dTm for the arm melting-temperature difference, or order for probe order.">?</span></label><select name="plot_color_by"><option value="gc">gc</option><option value="dtm">dTm</option><option value="order">order</option></select></div>
<div data-workflows="design manifest plot"><label>Plot DPI</label><input name="plot_dpi" type="number" min="72" value="300"></div>
<div data-workflows="design manifest plot"><label>Plot title</label><input name="plot_title" placeholder="Automatic title"></div>
<div data-workflows="plot"><label>Plot transcript length</label><input name="transcript_length" type="number" min="1" placeholder="Required only when it cannot be inferred"></div>
<div class="check" data-workflows="design manifest plot"><input name="no_probe_labels" type="checkbox"><label style="margin:0">Hide probe labels on maps</label><span class="help" tabindex="0" aria-label="Probe label help" data-help="Hide the pair numbers printed on maps while keeping the probe coordinates and data files unchanged.">?</span></div>
<div class="check" data-workflows="design manifest"><input name="no_oligo_qc" type="checkbox"><label style="margin:0">Skip final oligo QC</label><span class="help" tabindex="0" aria-label="Skip QC help" data-help="Skip the Primer3 check on final orderable oligos. This is not recommended for routine designs and cannot be combined with smart auto-curation.">?</span></div>
<div data-workflows="qc"><label class="field-label">QC workbook output path <span class="help" tabindex="0" aria-label="QC output help" data-help="Only used by the QC workflow. Leave blank to create the workbook automatically in the QC run folder.">?</span></label><input name="qc_output" placeholder="Automatic beside QC input"></div>
<div data-workflows="design manifest"><label class="field-label">QC stringency <span class="help" tabindex="0" aria-label="QC stringency help" data-help="Controls which QC reviews are allowed during smart selection: strict is most conservative, balanced is the default, and permissive allows more candidates during rescue.">?</span></label><select name="qc_stringency"><option value="balanced">balanced</option><option value="strict">strict</option><option value="permissive">permissive</option></select></div>
<div data-workflows="design manifest"><label class="field-label">Auto-curation plan <span class="help" tabindex="0" aria-label="Auto plan help" data-help="Controls how far smart fallback can relax tile, GC, Gibbs, and C/G-run settings. Standard is the usual choice; deep searches shorter tiles.">?</span></label><select name="auto_curate_plan"><option value="standard">standard</option><option value="conservative">conservative</option><option value="deep">deep</option></select></div>
<div data-workflows="design manifest"><label class="field-label">Target probe pairs <span class="help" tabindex="0" aria-label="Target probes help" data-help="Desired maximum number of non-overlapping pairs. It is a target, not a guarantee; sequence capacity and QC can limit the final count.">?</span></label><input name="target_probes" type="number" min="1" value="20"></div>
<div data-workflows="design manifest"><label class="field-label">Minimum acceptable pairs <span class="help" tabindex="0" aria-label="Minimum acceptable help" data-help="Smallest useful partial result used by smart curation when the requested target cannot be reached.">?</span></label><input name="min_acceptable_probes" type="number" min="1" value="12"></div>
<div data-workflows="design manifest"><label class="field-label">Maximum additional runs per phase <span class="help" tabindex="0" aria-label="Additional runs help" data-help="Limits how many distinct fallback design passes can run in each curation phase. More runs may improve yield but take longer.">?</span></label><input name="max_auto_runs" type="number" min="0" value="12"></div>
<div data-workflows="design manifest"><label class="field-label">Max probes per added tier <span class="help" tabindex="0" aria-label="Tier probes help" data-help="Candidate-reservoir size for each fallback tier. Leave blank to use the automatic value based on Target probe pairs.">?</span></label><input name="auto_curate_max_probes" type="number" min="1" placeholder="Dynamic default"></div>
<div data-workflows="design manifest"><label class="field-label">Coverage policy <span class="help" tabindex="0" aria-label="Coverage help" data-help="Balanced prefers even transcript-wide distribution after count and QC are tied. Choose qc-only when you want only the quality and technical objective.">?</span></label><select name="coverage_policy"><option value="balanced">balanced</option><option value="qc-only">qc-only</option></select></div>
</div></details></div>
<div class="full actions"><button type="submit">Run HCRProbeForge</button></div>
</div></form>
<div class="species-modal-backdrop hidden" id="species_modal" role="dialog" aria-modal="true" aria-labelledby="species_modal_title"><section class="species-modal-card"><button type="button" class="welcome-close" id="close_species_modal" aria-label="Close add species dialog">×</button><h2 id="species_modal_title">Add a new species</h2><p>Create a reusable species preset for an organism that is not in the built-in list. Required fields define transcript lookup and the user-facing name. Assembly details and a local genome FASTA are optional; you can supply them now or later in the Build a genome index workflow.</p><form id="species_form" enctype="multipart/form-data"><div class="grid"><div><label>Species alias <span class="required">*</span> <span class="help" tabindex="0" aria-label="Species alias help" data-help="A stable short identifier used by the CLI, index registry, cache, and output folders. Use letters, numbers, dots, underscores, or hyphens; for example octopus_v1.">?</span></label><input name="custom_key" required placeholder="e.g. octopus_v1"><div class="tip">This is not the scientific name; it is the internal alias.</div></div><div><label>Display name <span class="required">*</span> <span class="help" tabindex="0" aria-label="Display name help" data-help="The readable name shown in the webapp species menu, for example Octopus vulgaris.">?</span></label><input name="custom_display_name" required placeholder="e.g. Octopus vulgaris"></div><div class="full"><label>NCBI organism <span class="required">*</span> <span class="help" tabindex="0" aria-label="NCBI organism help" data-help="The exact organism text used in NCBI transcript searches. Use the taxon name recognized by NCBI, for example Octopus vulgaris.">?</span></label><input name="custom_organism" required placeholder="e.g. Octopus vulgaris"><div class="tip">This keeps transcript retrieval matched to the selected species.</div></div><div><label>Assembly accession <span class="optional">(optional)</span> <span class="help" tabindex="0" aria-label="Custom assembly accession help" data-help="A versioned NCBI GCF_ or GCA_ accession. If supplied, the Build index workflow can download and verify this assembly.">?</span></label><input name="custom_assembly_accession" placeholder="e.g. versioned GCF_ or GCA_ accession"></div><div><label>Assembly name <span class="optional">(optional)</span> <span class="help" tabindex="0" aria-label="Custom assembly name help" data-help="The expected assembly name returned by NCBI, or a descriptive name for your local genome FASTA.">?</span></label><input name="custom_assembly_name" placeholder="e.g. Octopus_v1"></div><div class="full"><label>Genome FASTA <span class="optional">(optional)</span> <span class="help" tabindex="0" aria-label="Custom genome FASTA help" data-help="Upload the genome FASTA you want the Build index workflow to use. HCRProbeForge stores a private copy in its reference-data directory; it does not alter the original file.">?</span></label><input name="custom_genome_fasta" type="file" accept=".fa,.fasta,.fna,.gz"><div class="tip">Use this when the genome is local or unavailable from NCBI. Transcript lookup still uses the NCBI organism above.</div></div></div><div class="species-modal-message error" id="species_modal_message"></div><div class="welcome-actions" style="margin-top:18px"><button type="button" class="secondary" id="cancel_species_modal">Cancel</button><button type="submit" id="save_species">Add species</button></div></form></section></div>
<div class="welcome-backdrop" id="welcome_modal" role="dialog" aria-modal="true" aria-labelledby="welcome_title"><section class="welcome-card"><button type="button" class="welcome-close" id="close_welcome" aria-label="Close getting started guide">×</button><h2 id="welcome_title">Getting started</h2><p>This app runs on your computer and saves the files in the results folder you choose. NCBI is contacted only when a workflow needs transcript or genome information.</p><ol><li><strong>Choose the species.</strong> Select the organism that matches your target. This keeps transcript lookup and genome specificity screening on the same reference.</li><li><strong>Check the index.</strong> “Index ready on this computer” means the selected Bowtie2 index is ready. Otherwise choose <em>Build a genome index</em>; a preset uses its compatible assembly automatically, while the optional accession and alias fields let you build another assembly.</li><li><strong>Choose your input.</strong> For one target, enter a gene and optionally an exact RefSeq accession, or provide one FASTA record. For many targets, upload a CSV, TSV, XLSX, XLSM, or text gene list. Choose one channel or <em>All channels</em> for B1–B5.</li><li><strong>Use the other workflows when needed.</strong> <em>Plot a probe table</em> creates a map from an existing probe table. <em>QC an oligo table</em> writes a QC workbook. Advanced settings expose the same thresholds and plot controls as the CLI.</li><li><strong>Find the results.</strong> The report shows the run folder and the key files. Selected pairs, order files, maps, and QC workbooks stay at the top of each target folder; intermediate files are kept in <code>details/</code>.</li></ol><div class="welcome-actions"><button type="button" class="secondary" id="welcome_done">Got it</button></div></section></div>
</main><script>
const mode=document.getElementById('mode');
const tabs=[...document.querySelectorAll('.workflow-tab')];
const show=(id,visible)=>{const node=document.getElementById(id);if(node){node.hidden=!visible;node.style.display=visible?'':'none';node.querySelectorAll('input,select,textarea').forEach(control=>{control.disabled=!visible;});}};
function alignAdvancedChecks(){
 const grid=document.getElementById('advanced_grid');
 if(!grid||window.matchMedia('(max-width:760px)').matches)return;
 const children=[...grid.children];
 children.filter(item=>item.classList.contains('check')).forEach(check=>{
   check.style.transform='';
   if(check.hidden)return;
   const checkTop=check.getBoundingClientRect().top;
   const target=children.find(item=>{
     if(item===check||item.hidden)return false;
     const rect=item.getBoundingClientRect();
     return Math.abs(rect.top-checkTop)<2 && item.querySelector('select,input:not([type="checkbox"]),textarea');
   });
   const control=target&&target.querySelector('select,input:not([type="checkbox"]),textarea');
   const checkbox=check.querySelector('input[type="checkbox"]');
   if(control&&checkbox){
     const controlRect=control.getBoundingClientRect();
     const checkboxRect=checkbox.getBoundingClientRect();
     const offset=controlRect.top+(controlRect.height-checkboxRect.height)/2-checkboxRect.top;
     check.style.transform=`translateY(${Math.round(offset)}px)`;
   }
 });
}
function updateMode(){
 const value=mode.value;
 tabs.forEach(tab=>{const active=tab.dataset.mode===value;tab.classList.toggle('active',active);tab.setAttribute('aria-selected',active?'true':'false');});
 show('reference_input',true);show('channel_input',value==='design'||value==='manifest');
 show('design_inputs',value==='design');show('accession_input',value==='design');show('fasta_input',value==='design');show('target_type_input',value==='design');
 show('manifest_input',value==='manifest');show('existing_input',value==='plot'||value==='qc');show('index_input',value==='index');show('index_local_fasta_input',value==='index');
 show('output_input',value!=='index');show('folder_input',value!=='index');show('auto_input',value==='design'||value==='manifest');show('design_note',value==='design'||value==='manifest');
 const advanced=document.getElementById('advanced_settings');
 const advancedGrid=document.getElementById('advanced_grid');
 const isDesign=value==='design'||value==='manifest';
 advanced.hidden=value==='index';
 advanced.style.display=value==='index'?'none':'';
 advancedGrid.dataset.mode=value;
 [...advancedGrid.children].filter(item=>item.dataset.workflows).forEach(item=>{
   const visible=(item.dataset.workflows||'').split(/\s+/).includes(value);
   item.hidden=!visible;
   item.setAttribute('aria-hidden',visible?'false':'true');
   item.querySelectorAll('[name]').forEach(control=>{control.disabled=!visible;});
 });
 document.getElementById('auto').disabled=!isDesign;
 requestAnimationFrame(alignAdvancedChecks);
}
function updateTargetControls(){
 const target=document.getElementById('target_type');
 const details=document.getElementById('premrna_details');
 if(!target||!details)return;
 const preOptions=target.querySelectorAll('option[value^="pre-mrna"]');
 const currentSpecies=document.getElementById('species');
const selected=currentSpecies&&currentSpecies.options[currentSpecies.selectedIndex];
const hasPreset=!!(selected&&selected.value&&selected.value!=='__add_new_species__');
preOptions.forEach(option=>{option.hidden=!hasPreset;});
if(!hasPreset&&target.value.startsWith('pre-mrna'))target.value='mature';
const isPre=target.value.startsWith('pre-mrna')&&hasPreset&&mode.value==='design';
const isIntronic=target.value==='pre-mrna'&&isPre;
const targetTip=document.getElementById('target_type_tip');
if(targetTip){
  targetTip.textContent=target.value==='pre-mrna'
    ? 'Pre-mRNA · intronic regions: use the full genomic transcript, but keep every probe pair inside the selected introns. Leave the field blank for all annotated introns, or enter 1,2 for specific introns.'
    : target.value==='pre-mrna-whole'
      ? 'Pre-mRNA · whole genomic transcript: scan the full genomic transcript. Probes may fall in exons, introns, or across exon–intron boundaries.'
      : 'Mature transcript: design on the processed, spliced RNA. Choose a Pre-mRNA mode when you need unspliced nuclear RNA; it requires a matching genomic annotation database.';
}
details.hidden=!isIntronic;
details.querySelectorAll('input').forEach(control=>{control.disabled=!isIntronic;});
const annotationStatus=document.getElementById('premrna_annotation_status');
if(annotationStatus){const ready=selected&&selected.dataset.annotationReady==='true';annotationStatus.hidden=!isPre;annotationStatus.textContent=ready?'Pre-mRNA database ready':'Pre-mRNA database will be built on first Pre-mRNA run';annotationStatus.classList.toggle('ready',!!ready);annotationStatus.classList.toggle('pending',!ready);}
const auto=document.getElementById('auto');
if(auto&&mode.value==='design'){auto.disabled=false;auto.title=isIntronic?'Optional: if enabled, curation distributes pairs across eligible introns.':isPre?'Optional: if enabled, curation searches fallback tiers across the full genomic target.':'';}
}
tabs.forEach(tab=>tab.addEventListener('click',()=>{mode.value=tab.dataset.mode;updateMode();updateTargetControls();}));mode.addEventListener('change',()=>{updateMode();updateTargetControls();});window.addEventListener('resize',()=>requestAnimationFrame(alignAdvancedChecks));updateMode();updateTargetControls();
const species=document.getElementById('species');const organism=document.getElementById('organism');const indexStatus=document.getElementById('index_status');const indexPresetReference=document.getElementById('index_preset_reference');const addSpeciesValue='__add_new_species__';let selectedSpecies=species.value;function updateOrganism(){const option=species.options[species.selectedIndex];organism.value=option?option.dataset.organism||'':'';if(indexStatus){const ready=option&&option.dataset.indexReady==='true';indexStatus.textContent=ready?'Index ready on this computer':'Index not registered yet';indexStatus.classList.toggle('ready',!!ready);}if(indexPresetReference){const assembly=option&&option.dataset.assemblyName;const accession=option&&option.dataset.assemblyAccession;indexPresetReference.textContent=assembly&&accession?`The selected preset is ${assembly} (${accession}).`:'The selected preset supplies its compatible assembly.';}updateTargetControls();}species.addEventListener('change',()=>{if(species.value===addSpeciesValue){species.value=selectedSpecies;updateOrganism();openSpeciesModal();return;}selectedSpecies=species.value;updateOrganism();});updateOrganism();
const welcome=document.getElementById('welcome_modal');const closeWelcome=()=>{welcome.classList.add('hidden');try{window.localStorage.setItem('hcrprobeforge_welcome_seen','1');}catch(e){}};document.getElementById('close_welcome').addEventListener('click',closeWelcome);document.getElementById('welcome_done').addEventListener('click',closeWelcome);document.getElementById('open_welcome').addEventListener('click',()=>{welcome.classList.remove('hidden');});welcome.addEventListener('click',event=>{if(event.target===welcome)closeWelcome();});try{if(window.localStorage.getItem('hcrprobeforge_welcome_seen')==='1')welcome.classList.add('hidden');}catch(e){}
document.getElementById('inspect').addEventListener('click',async()=>{const box=document.getElementById('inspection');box.style.display='block';box.textContent='Inspecting input…';try{const r=await fetch('/inspect',{method:'POST',body:new FormData(document.querySelector('form'))});const data=await r.json();if(!r.ok)throw new Error(data.error||'Inspection failed');box.textContent=data.message;}catch(e){box.style.background='var(--rose)';box.textContent=e.message;}});
const speciesModal=document.getElementById('species_modal');const speciesForm=document.getElementById('species_form');const speciesMessage=document.getElementById('species_modal_message');const openSpeciesModal=()=>{speciesMessage.style.display='none';speciesMessage.textContent='';speciesModal.classList.remove('hidden');speciesForm.elements.custom_key.focus();};const closeSpeciesModal=()=>{speciesModal.classList.add('hidden');};document.getElementById('close_species_modal').addEventListener('click',closeSpeciesModal);document.getElementById('cancel_species_modal').addEventListener('click',closeSpeciesModal);speciesModal.addEventListener('click',event=>{if(event.target===speciesModal)closeSpeciesModal();});speciesForm.addEventListener('submit',async event=>{event.preventDefault();const button=document.getElementById('save_species');button.disabled=true;button.textContent='Adding…';speciesMessage.style.display='none';try{const response=await fetch('/register-species',{method:'POST',body:new FormData(speciesForm)});const data=await response.json();if(!response.ok)throw new Error(data.error||'Could not add species');const option=document.createElement('option');option.value=data.key;option.textContent=data.display_name;option.dataset.organism=data.scientific_name;option.dataset.assemblyName=data.assembly_name||'';option.dataset.assemblyAccession=data.assembly_accession||'';option.dataset.indexReady='false';const addOption=species.querySelector('option[value="'+addSpeciesValue+'"]');species.insertBefore(option,addOption);species.value=data.key;selectedSpecies=data.key;updateOrganism();closeSpeciesModal();speciesForm.reset();}catch(error){speciesMessage.style.display='block';speciesMessage.textContent=error.message;}finally{button.disabled=false;button.textContent='Add species';}});
</script></body></html>
"""


def _parse_multipart(handler: BaseHTTPRequestHandler) -> tuple[dict[str, str], list[tuple[str, str, bytes]]]:
    length = int(handler.headers.get("Content-Length", "0"))
    content_type = handler.headers.get("Content-Type", "")
    raw = handler.rfile.read(length)
    message = BytesParser(policy=default).parsebytes(
        f"Content-Type: {content_type}\r\nMIME-Version: 1.0\r\n\r\n".encode() + raw
    )
    fields: dict[str, str] = {}
    uploads: list[tuple[str, str, bytes]] = []
    if not message.is_multipart():
        return fields, uploads
    for part in message.iter_parts():
        name = part.get_param("name", header="content-disposition") or ""
        filename = part.get_filename()
        payload = part.get_payload(decode=True) or b""
        if filename and payload:
            uploads.append((name, Path(filename).name, payload))
        elif name:
            fields[name] = payload.decode(part.get_content_charset() or "utf-8", errors="replace")
    return fields, uploads


def _upload_or_path(
    fields: dict[str, str],
    uploads: list[tuple[str, str, bytes]],
    field: str,
    upload_field: str,
    job_dir: Path,
) -> Path | None:
    for name, filename, payload in uploads:
        if name == upload_field:
            target = job_dir / (core.safe_name(Path(filename).name) or "input_file")
            target.write_bytes(payload)
            return target
    value = fields.get(field, "").strip()
    if value:
        path = Path(value).expanduser().resolve()
        if not path.is_file():
            raise ValueError(f"file not found: {path}")
        return path
    return None


def _preserve_uploaded_input(
    path: Path | None,
    uploads: list[tuple[str, str, bytes]],
    upload_field: str,
    persistent_root: Path,
) -> Path | None:
    """Keep an uploaded input stable after its request staging directory closes.

    User-supplied local paths are returned unchanged. Uploaded inputs are copied
    below the workflow output so run metadata and index metadata never point at
    a deleted temporary file.
    """
    if path is None or not any(name == upload_field for name, _, _ in uploads):
        return path
    input_root = persistent_root / ".hcrprobeforge-inputs" / uuid4().hex
    input_root.mkdir(parents=True, exist_ok=True)
    preserved = input_root / (core.safe_name(path.name) or "input_file")
    shutil.copyfile(path, preserved)
    return preserved


def _add_value(argv: list[str], fields: dict[str, str], field: str, option: str) -> None:
    value = fields.get(field, "").strip()
    if value:
        argv.extend([option, value])


def _add_flag(argv: list[str], fields: dict[str, str], field: str, option: str) -> None:
    if field in fields:
        argv.append(option)


def _common_design_args(fields: dict[str, str]) -> list[str]:
    argv: list[str] = []
    for field, option in (
        ("organism", "--organism"), ("species", "--species"), ("index", "--index"), ("email", "--email"),
        ("api_key", "--api-key"), ("tile_size", "--tile-size"), ("min_gc", "--min-gc"), ("max_gc", "--max-gc"),
        ("min_gibbs", "--min-gibbs"), ("max_gibbs", "--max-gibbs"), ("target_gibbs", "--target-gibbs"),
        ("max_run_mismatches", "--max-run-mismatches"), ("max_probes", "--max-probes"),
        ("num_hits_allowed", "--num-hits-allowed"), ("dtm_max", "--dtm-max"),
        ("target_probes", "--target-probes"), ("min_acceptable_probes", "--min-acceptable-probes"),
        ("max_auto_runs", "--max-auto-runs"), ("auto_curate_max_probes", "--auto-curate-max-probes"),
        ("target_type", "--target-type"), ("premrna_introns", "--premrna-introns"),
    ):
        _add_value(argv, fields, field, option)
    for field, option in (
        ("all_transcripts", "--all-transcripts"), ("dtm_filter", "--dtm-filter"), ("no_genomemask", "--no-genomemask"),
        ("auto_curate", "--auto-curate-if-needed"), ("no_oligo_qc", "--no-oligo-qc"),
    ):
        _add_flag(argv, fields, field, option)
    for field, option in (
        ("transcript_policy", "--transcript-policy"), ("qc_stringency", "--qc-stringency"),
        ("auto_curate_plan", "--auto-curate-plan"), ("coverage_policy", "--coverage-policy"),
    ):
        _add_value(argv, fields, field, option)
    argv += _plot_render_args(fields, include_transcript_length=False)
    return argv


def _plot_render_args(fields: dict[str, str], *, include_transcript_length: bool) -> list[str]:
    argv: list[str] = []
    plot_fields = [
        ("plot_theme", "--plot-theme"),
        ("plot_color_by", "--plot-color-by"),
        ("plot_dpi", "--plot-dpi"),
        ("plot_title", "--plot-title"),
    ]
    if include_transcript_length:
        plot_fields.append(("transcript_length", "--transcript-length"))
    for field, option in plot_fields:
        _add_value(argv, fields, field, option)
    _add_flag(argv, fields, "no_probe_labels", "--no-probe-labels")
    return argv


def _plot_args(fields: dict[str, str]) -> list[str]:
    return _plot_render_args(fields, include_transcript_length=True)


def _output_path(fields: dict[str, str]) -> Path:
    return Path(fields.get("outdir", "hcr_results").strip() or "hcr_results").expanduser().resolve()


def _register_species_submission(
    fields: dict[str, str],
    uploads: list[tuple[str, str, bytes]],
) -> dict[str, object]:
    """Validate and persist one species dialog submission."""
    with tempfile.TemporaryDirectory(prefix="hcrprobeforge-species-") as temporary:
        return _register_species_submission_in_directory(fields, uploads, Path(temporary))


def _register_species_submission_in_directory(
    fields: dict[str, str],
    uploads: list[tuple[str, str, bytes]],
    job_dir: Path,
) -> dict[str, object]:
    genome = _upload_or_path(fields, uploads, "custom_genome_fasta_path", "custom_genome_fasta", job_dir)
    preset = references.register_custom_species(
        fields.get("custom_key", ""),
        fields.get("custom_display_name", ""),
        fields.get("custom_organism", ""),
        assembly_accession=fields.get("custom_assembly_accession", "") or None,
        assembly_name=fields.get("custom_assembly_name", "") or None,
        genome_fasta=genome,
    )
    return {
        "key": preset.key,
        "display_name": preset.display_name,
        "scientific_name": preset.scientific_name,
        "assembly_name": preset.assembly_name,
        "assembly_accession": preset.assembly_accession,
        "local_genome_fasta": bool(preset.local_genome_fasta),
    }


def _delete_species_submission(fields: dict[str, str]) -> dict[str, str]:
    key = fields.get("species", "").strip()
    plan = references.species_deletion_plan(key)
    references.delete_species(key)
    return {
        "key": references.canonical_species(key),
        "removed": "true",
        "builtin": "true" if plan.get("builtin") else "false",
        "preset_protected": "true" if plan.get("preset_protected") else "false",
        "deleted_index_files": str(len(plan.get("index_files", []))),
        "preserved_shared_index_files": str(len(plan.get("preserved_index_files", []))),
    }


def _invoke_with_cache(
    argv: list[str],
    *,
    output_base: Path,
    species: str,
    workflow: str,
    progress_callback: ProgressCallback | None = None,
) -> tuple[int, str, str]:
    """Invoke the core CLI with the webapp's project-level cache selected."""
    cache_root = references.species_cache_root(output_base, species, workflow)
    # The core creates cache subdirectories only when it has data to cache.
    original_cache = os.environ.get("HCRPROBEFORGE_CACHE_DIR")
    original_species = os.environ.get("HCRPROBEFORGE_SPECIES")
    original_workflow = os.environ.get("HCRPROBEFORGE_WORKFLOW")
    os.environ["HCRPROBEFORGE_CACHE_DIR"] = str(cache_root)
    os.environ["HCRPROBEFORGE_SPECIES"] = references.get_species_preset(species).key if references.get_species_preset(species) else species
    os.environ["HCRPROBEFORGE_WORKFLOW"] = workflow
    try:
        return batch._invoke_core(
            argv,
            capture=True,
            progress_callback=progress_callback,
        )
    finally:
        if original_cache is None:
            os.environ.pop("HCRPROBEFORGE_CACHE_DIR", None)
        else:
            os.environ["HCRPROBEFORGE_CACHE_DIR"] = original_cache
        if original_species is None:
            os.environ.pop("HCRPROBEFORGE_SPECIES", None)
        else:
            os.environ["HCRPROBEFORGE_SPECIES"] = original_species
        if original_workflow is None:
            os.environ.pop("HCRPROBEFORGE_WORKFLOW", None)
        else:
            os.environ["HCRPROBEFORGE_WORKFLOW"] = original_workflow


def _species_options_html(
    *,
    index_ready: dict[str, bool] | None = None,
    annotation_ready: dict[str, bool] | None = None,
    installed_references: list[dict[str, object]] | None = None,
    reference_metadata: dict[str, dict[str, object]] | None = None,
    annotation_structural_check: bool = True,
) -> str:
    """Render user presets and registered alternate assemblies."""
    custom_presets = tuple(
        preset
        for preset in references.custom_species_presets()
        if references.canonical_species(preset.key) not in references.REMOVED_BUILTIN_ALIASES
    )
    preset_keys = {preset.key for preset in references.supported_species()}
    custom_keys = {preset.key for preset in custom_presets}
    custom_display_names = {references.canonical_species(preset.display_name) for preset in custom_presets}
    custom_organisms = {references.canonical_species(preset.scientific_name) for preset in custom_presets}
    builtin_reference_identity = {
        (
            references.canonical_species(preset.scientific_name),
            references.canonical_species(preset.assembly_name),
            str(preset.assembly_accession or "").strip().upper(),
        )
        for preset in references.SPECIES_PRESETS.values()
    }
    builtin_reference_pairs = {
        (
            references.canonical_species(preset.scientific_name),
            references.canonical_species(preset.assembly_name),
        )
        for preset in references.SPECIES_PRESETS.values()
    }
    options: list[str] = []
    seen: set[str] = set()
    for preset in custom_presets:
        seen.add(preset.key)
        preset_metadata = (
            reference_metadata.get(references.canonical_species(preset.key))
            if reference_metadata is not None
            else None
        )
        preset_index_ready = (
            index_ready.get(preset.key)
            if index_ready is not None and preset.key in index_ready
            else references.registered_index_is_ready(preset.key, metadata=preset_metadata)
        )
        preset_annotation_ready = (
            annotation_ready.get(preset.key)
            if annotation_ready is not None and preset.key in annotation_ready
            else _annotation_database_ready(
                preset.key,
                metadata=preset_metadata,
                structural_check=annotation_structural_check,
            )
        )
        options.append(
            f'<option data-custom="true" data-local-fasta="{str(bool(preset.local_genome_fasta)).lower()}" data-organism="{html.escape(preset.scientific_name, quote=True)}" '
            f'data-assembly-name="{html.escape(preset.assembly_name, quote=True)}" '
            f'data-assembly-accession="{html.escape(preset.assembly_accession, quote=True)}" '
            f'data-index-ready="{str(bool(preset_index_ready)).lower()}" '
            f'data-annotation-ready="{str(bool(preset_annotation_ready)).lower()}" '
            f'value="{html.escape(preset.key, quote=True)}">{html.escape(preset.display_name)}</option>'
        )
    for row in installed_references if installed_references is not None else references.list_installed_references():
        if row.get("status") != "ready":
            continue
        alias = core.safe_name(str(row.get("index_species_alias") or row.get("species") or ""))
        if not alias or alias in preset_keys or alias in seen:
            continue
        if references.canonical_species(alias) in references.REMOVED_BUILTIN_ALIASES:
            # Do not expose the removed legacy Xenopus alias from stale
            # HCRProbeDesign registrations or old reference metadata.
            continue
        scientific = str(row.get("scientific_name") or "").strip()
        row_species = references.canonical_species(str(row.get("species") or ""))
        row_display = references.canonical_species(str(row.get("display_name") or ""))
        row_organism = references.canonical_species(scientific)
        row_assembly = references.canonical_species(
            str(row.get("assembly") or row.get("assembly_name") or "")
        )
        row_accession = str(
            row.get("assembly_accession") or row.get("accession") or ""
        ).strip().upper()
        is_duplicate_builtin_reference = (
            row_organism,
            row_assembly,
            row_accession,
        ) in builtin_reference_identity
        is_duplicate_builtin_reference = is_duplicate_builtin_reference or (
            bool(row_organism)
            and bool(row_assembly)
            and (row_organism, row_assembly) in builtin_reference_pairs
        )
        if (
            row_species in custom_keys
            or row_display in custom_display_names
            or row_organism in custom_organisms
            or is_duplicate_builtin_reference
        ):
            # A built-in or custom preset is already represented by its
            # user-facing option. Do not leak an alternate internal index
            # registration beside it.
            continue
        seen.add(alias)
        assembly = str(row.get("assembly") or "").strip()
        accession = str(row.get("assembly_accession") or "").strip()
        detail = " · ".join(part for part in (assembly, accession) if part)
        label = " · ".join(part for part in (alias, scientific, detail) if part)
        installed_index_ready = (
            index_ready.get(alias)
            if index_ready is not None and alias in index_ready
            else references.registered_index_is_ready(alias)
        )
        installed_annotation_ready = (
            annotation_ready.get(alias)
            if annotation_ready is not None and alias in annotation_ready
            else _annotation_database_ready(
                alias,
                metadata=row,
                structural_check=annotation_structural_check,
            )
        )
        options.append(
            f'<option data-organism="{html.escape(scientific, quote=True)}" '
            f'data-index-ready="{str(bool(installed_index_ready)).lower()}" '
            f'data-annotation-ready="{str(bool(installed_annotation_ready)).lower()}" '
            f'value="{html.escape(alias, quote=True)}">{html.escape(label)}</option>'
        )
    return "".join(options) + '<option value="__add_new_species__">＋ Add a new species…</option>'


def _annotation_database_ready(
    species: str,
    *,
    metadata: dict[str, object] | None = None,
    structural_check: bool = True,
) -> bool:
    """Return whether a matching Pre-mRNA lookup database is available.

    The initial setup page uses ``structural_check=False`` so it returns
    promptly after an index build. A completed metadata record plus the final
    database path is authoritative there because database creation is atomic.
    The selected design path and the live status endpoint use the full
    read-only SQLite schema check before consuming a database.
    """
    try:
        from . import premrna

        preset = references.get_species_preset(species)
        metadata = metadata if metadata is not None else (references.find_installed_reference(species) or {})
        # Installed alternate assemblies (for example c_elegans) are not
        # built-in SpeciesPreset objects. Their reference.json is still the
        # authoritative source for the index and annotation paths.
        species_key = preset.key if preset is not None else str(
            metadata.get("index_species_alias") or metadata.get("species") or species
        )
        if preset is None and not metadata:
            return False
        if not references.registered_index_is_ready(species_key, metadata=metadata or None):
            return False
        candidates: list[Path] = []
        configured = str(metadata.get("annotation_database_path") or "").strip()
        if configured:
            candidates.append(Path(configured).expanduser())
        reference_dir = str(metadata.get("reference_data_directory") or "").strip()
        if reference_dir:
            candidates.append(Path(reference_dir).expanduser() / "annotation.sqlite")
        if preset is not None:
            assembly = core.safe_name(preset.assembly_name or preset.assembly_accession or "assembly") or "assembly"
            candidates.append(references.species_data_root() / preset.key / assembly / "annotation.sqlite")
        if not structural_check:
            status = str(metadata.get("annotation_database_status") or "").casefold()
            if metadata and status not in {"", "ready"}:
                return False
            return any(path.is_file() for path in candidates)
        return any(premrna.annotation_database_is_structurally_ready(path) for path in candidates)
    except (OSError, TypeError, ValueError, sqlite3.Error):
        return False


WEBAPP_ENHANCEMENTS = r"""
const removeSpeciesButton=document.createElement('button');removeSpeciesButton.type='button';removeSpeciesButton.id='remove_species';removeSpeciesButton.className='secondary species-remove';removeSpeciesButton.hidden=true;removeSpeciesButton.textContent='Remove this custom preset';document.getElementById('species_input').appendChild(removeSpeciesButton);
const speciesPicker=document.createElement('div');speciesPicker.className='species-picker';species.parentNode.insertBefore(speciesPicker,species);speciesPicker.appendChild(species);speciesPicker.appendChild(removeSpeciesButton);removeSpeciesButton.textContent='×';removeSpeciesButton.setAttribute('aria-label','Remove selected custom species preset');removeSpeciesButton.title='Remove selected custom species preset';removeSpeciesButton.dataset.help='Remove this custom species preset. Built-in presets are protected; custom presets can be removed whether or not an index exists.';
const floatingHelp=document.createElement('div');floatingHelp.className='floating-help';floatingHelp.setAttribute('role','tooltip');floatingHelp.hidden=true;document.body.appendChild(floatingHelp);let floatingHelpOwner=null;function hideFloatingHelp(){floatingHelp.hidden=true;floatingHelpOwner=null;}function showFloatingHelp(node){const text=node.dataset.help||'';if(!text)return;floatingHelp.textContent=text;floatingHelp.hidden=false;floatingHelpOwner=node;const rect=node.getBoundingClientRect();const width=Math.min(360,window.innerWidth-24);const left=Math.max(12,Math.min(rect.left,window.innerWidth-width-12));floatingHelp.style.width=width+'px';const height=floatingHelp.offsetHeight;const top=rect.bottom+8+height<=window.innerHeight-12?rect.bottom+8:Math.max(12,rect.top-height-8);floatingHelp.style.left=left+'px';floatingHelp.style.top=top+'px';}document.querySelectorAll('.species-modal-card .help').forEach(node=>{node.addEventListener('mouseenter',()=>showFloatingHelp(node));node.addEventListener('mouseleave',hideFloatingHelp);node.addEventListener('focus',()=>showFloatingHelp(node));node.addEventListener('blur',hideFloatingHelp);});window.addEventListener('resize',()=>{if(floatingHelpOwner)showFloatingHelp(floatingHelpOwner);});
function syncCustomSpeciesControls(){const option=species.options[species.selectedIndex];const ready=option&&option.dataset.indexReady==='true';const custom=option&&option.dataset.custom==='true';removeSpeciesButton.hidden=!custom;removeSpeciesButton.dataset.help=custom?'Remove this custom species preset, its saved HCRProbeForge data, and its unshared Bowtie2 index files. Built-in presets are protected.':'Remove this custom species preset.';const assembly=option&&option.dataset.assemblyName;const accession=option&&option.dataset.assemblyAccession;const local=option&&option.dataset.localFasta==='true';const indexNote=document.querySelector('#index_input .note');if(indexNote&&indexNote.childNodes[1])indexNote.childNodes[1].textContent=custom?' select a species above. The NCBI organism identifies transcripts, but a custom preset still needs a versioned assembly accession or local genome FASTA.':' select a species above and the preset will supply the current compatible NCBI assembly. ';if(indexPresetReference&&!ready){indexPresetReference.textContent=assembly&&accession?`The selected custom preset is ${assembly} (${accession}).`:local?'The selected custom preset has a saved genome FASTA that will be used automatically.':'This custom preset has no genome source yet. Enter an NCBI assembly accession or upload a genome FASTA below.';}}
const deleteSpeciesModal=document.getElementById('delete_species_modal');const deleteSpeciesText=document.getElementById('delete_species_text');const deleteSpeciesAssets=document.getElementById('delete_species_assets');const deleteSpeciesMessage=document.getElementById('delete_species_message');const deleteSpeciesCheck=document.getElementById('confirm_delete_species');const deleteSpeciesConfirm=document.getElementById('confirm_delete_species_button');let pendingSpeciesOption=null;const closeDeleteSpecies=()=>{deleteSpeciesModal.classList.add('hidden');pendingSpeciesOption=null;deleteSpeciesCheck.checked=false;deleteSpeciesCheck.disabled=false;deleteSpeciesConfirm.disabled=true;deleteSpeciesMessage.style.display='none';deleteSpeciesMessage.textContent='';deleteSpeciesAssets.textContent='';};const openDeleteSpecies=async option=>{pendingSpeciesOption=option;const requestedKey=option.value;deleteSpeciesText.textContent=`Reviewing what will be removed for “${option.textContent}”…`;deleteSpeciesAssets.textContent='Loading the exact private-data, config, and Bowtie2-file cleanup list…';deleteSpeciesCheck.checked=false;deleteSpeciesCheck.disabled=true;deleteSpeciesConfirm.disabled=true;deleteSpeciesMessage.style.display='none';deleteSpeciesModal.classList.remove('hidden');try{const response=await fetch('/delete-species-preview?species='+encodeURIComponent(requestedKey),{cache:'no-store'});const plan=await response.json();if(!response.ok)throw new Error(plan.error||'Could not inspect species cleanup');if(pendingSpeciesOption!==option)return;const fasta=plan.saved_fasta?`Saved FASTA: ${plan.saved_fasta}. `:'No saved private FASTA. ';const metadata=(plan.reference_metadata||[]).length?`HCRProbeForge metadata (${plan.reference_metadata.length}): ${plan.reference_metadata.join(', ')}. `:'No HCRProbeForge reference metadata. ';const config=(plan.config_files||[]).length?`HCRProbeDesign config entries (${plan.config_aliases_to_remove.length}) will be removed from: ${plan.config_files.join(', ')}. `:'No stale HCRProbeDesign config entry was found. ';const files=plan.index_files||[];const indexes=files.length?`Bowtie2 index/temp files (${files.length}) that will be deleted: ${files.join(', ')}.`:`No matching Bowtie2 index files were found under ${plan.index_root}.`;const folders=(plan.index_directories||[]).length?` Empty index folders that will be removed: ${plan.index_directories.join(', ')}.`:'';const shared=(plan.preserved_index_files||[]).length?` Shared index files retained for safety (${plan.preserved_index_files.length}): ${plan.preserved_index_files.join(', ')}.`:'';const unsafe=(plan.unsafe_index_prefixes||[]).length?` Index prefixes outside the managed HCRProbeDesign folder are retained: ${plan.unsafe_index_prefixes.join(', ')}.`:'';deleteSpeciesText.textContent=`This removes “${option.textContent}” from the species list, its saved preset data, stale config registrations, and the listed unshared index files. Built-in presets are protected.`;deleteSpeciesAssets.textContent=fasta+metadata+config+indexes+folders+shared+unsafe;deleteSpeciesCheck.disabled=false;deleteSpeciesCheck.focus();}catch(error){if(pendingSpeciesOption===option){deleteSpeciesText.textContent='The cleanup list could not be inspected, so deletion is disabled.';deleteSpeciesAssets.textContent='';deleteSpeciesMessage.style.display='block';deleteSpeciesMessage.textContent=error.message;}}};document.getElementById('close_delete_species').addEventListener('click',closeDeleteSpecies);document.getElementById('cancel_delete_species').addEventListener('click',closeDeleteSpecies);deleteSpeciesModal.addEventListener('click',event=>{if(event.target===deleteSpeciesModal)closeDeleteSpecies();});deleteSpeciesCheck.addEventListener('change',()=>{deleteSpeciesConfirm.disabled=!deleteSpeciesCheck.checked;});deleteSpeciesConfirm.addEventListener('click',async()=>{const option=pendingSpeciesOption;if(!option||!deleteSpeciesCheck.checked)return;const body=new FormData();body.append('species',option.value);deleteSpeciesConfirm.disabled=true;deleteSpeciesCheck.disabled=true;try{const response=await fetch('/delete-species',{method:'POST',body});const data=await response.json();if(!response.ok)throw new Error(data.error||'Could not remove species');option.remove();species.selectedIndex=0;selectedSpecies=species.value;updateOrganism();syncCustomSpeciesControls();closeDeleteSpecies();}catch(error){deleteSpeciesMessage.style.display='block';deleteSpeciesMessage.textContent=error.message;deleteSpeciesConfirm.disabled=false;}finally{deleteSpeciesCheck.disabled=false;}});species.addEventListener('change',()=>setTimeout(syncCustomSpeciesControls,0));new MutationObserver(syncCustomSpeciesControls).observe(species,{childList:true});syncCustomSpeciesControls();removeSpeciesButton.addEventListener('click',()=>{const option=species.options[species.selectedIndex];if(option&&option.dataset.custom==='true')openDeleteSpecies(option);});
['mouseenter','focus'].forEach(eventName=>removeSpeciesButton.addEventListener(eventName,()=>showFloatingHelp(removeSpeciesButton)));['mouseleave','blur'].forEach(eventName=>removeSpeciesButton.addEventListener(eventName,hideFloatingHelp));
const indexGrid=document.querySelector('#index_input .grid');const rebuildControl=indexGrid&&indexGrid.querySelector('[name="index_force"]')?.closest('.check');if(rebuildControl){const openIndexFolderButton=document.createElement('button');openIndexFolderButton.type='button';openIndexFolderButton.className='secondary index-folder-button';openIndexFolderButton.textContent='Open index folder';openIndexFolderButton.setAttribute('aria-label','Open the folder where genome indices are built');openIndexFolderButton.dataset.help='Open the HCRProbeDesign data folder that contains the built genome indices.';rebuildControl.appendChild(openIndexFolderButton);['mouseenter','focus'].forEach(eventName=>openIndexFolderButton.addEventListener(eventName,()=>showFloatingHelp(openIndexFolderButton)));['mouseleave','blur'].forEach(eventName=>openIndexFolderButton.addEventListener(eventName,hideFloatingHelp));openIndexFolderButton.addEventListener('click',async()=>{openIndexFolderButton.disabled=true;try{const response=await fetch('/open-index-folder',{method:'POST'});const data=await response.json();if(!response.ok)throw new Error(data.error||'Could not open the index folder');}catch(error){window.alert(error.message);}finally{openIndexFolderButton.disabled=false;}});}
async function refreshIndexStatus(){const option=species.options[species.selectedIndex];if(!option||!option.value||option.value===addSpeciesValue)return;try{const response=await fetch('/index-status?species='+encodeURIComponent(option.value),{cache:'no-store'});const data=await response.json();if(!response.ok)throw new Error(data.error||'Index status unavailable');option.dataset.indexReady=data.ready?'true':'false';option.dataset.annotationReady=data.annotation_ready?'true':'false';updateOrganism();syncCustomSpeciesControls();}catch(error){/* A transient status failure must not disrupt an active form. */}}
refreshIndexStatus();setInterval(refreshIndexStatus,3000);window.addEventListener('focus',refreshIndexStatus);window.addEventListener('pageshow',()=>refreshIndexStatus());document.addEventListener('visibilitychange',()=>{if(!document.hidden)refreshIndexStatus();});
const runForm=document.querySelector('form[action="/run"]');if(runForm)runForm.addEventListener('submit',()=>{window.hcrNavigating=true;});
const indexAnnotationModal=document.getElementById('index_annotation_modal');const indexAnnotationChoice=document.getElementById('index_annotation_database_choice');const indexAnnotationSpecies=document.getElementById('index_annotation_species');const closeIndexAnnotation=()=>{indexAnnotationModal.classList.add('hidden');indexAnnotationChoice.checked=false;window.hcrNavigating=false;};const openIndexAnnotation=()=>{const option=species.options[species.selectedIndex];const name=option?.textContent?.trim()||'the selected species';const assembly=option?.dataset?.assemblyName||'';const accession=option?.dataset?.assemblyAccession||'';indexAnnotationSpecies.textContent=`The selected reference is ${name}${assembly?` · ${assembly}`:''}${accession?` (${accession})`:''}.`;indexAnnotationChoice.checked=false;indexAnnotationModal.classList.remove('hidden');indexAnnotationChoice.focus();};document.getElementById('close_index_annotation').addEventListener('click',closeIndexAnnotation);document.getElementById('cancel_index_annotation').addEventListener('click',closeIndexAnnotation);indexAnnotationModal.addEventListener('click',event=>{if(event.target===indexAnnotationModal)closeIndexAnnotation();});if(runForm)runForm.addEventListener('submit',event=>{if(mode.value!=='index'||runForm.dataset.indexSubmitArmed==='1')return;event.preventDefault();window.hcrNavigating=false;openIndexAnnotation();});document.getElementById('continue_index_build').addEventListener('click',()=>{let hidden=runForm.querySelector('input[name="index_annotation_database"]');if(!hidden){hidden=document.createElement('input');hidden.type='hidden';hidden.name='index_annotation_database';runForm.appendChild(hidden);}hidden.value=indexAnnotationChoice.checked?'1':'0';runForm.dataset.indexSubmitArmed='1';window.hcrNavigating=true;indexAnnotationModal.classList.add('hidden');runForm.submit();});document.querySelectorAll('a[href^="/"]').forEach(link=>link.addEventListener('click',()=>{window.hcrNavigating=true;}));
"""


BROWSER_LIFECYCLE_SCRIPT = r'''<script>
window.hcrNavigating=false;
document.querySelectorAll('a[href^="/"]').forEach(link=>link.addEventListener('click',()=>{window.hcrNavigating=true;}));
document.querySelectorAll('a[href^="/?new_run="]').forEach(link=>link.addEventListener('click',()=>{link.textContent='Opening setup…';link.setAttribute('aria-busy','true');link.style.pointerEvents='none';}));
document.querySelectorAll('form').forEach(form=>form.addEventListener('submit',()=>{window.hcrNavigating=true;}));
</script>'''


def _browser_lifecycle_script(request_token: str) -> str:
    """Return navigation bookkeeping without tying tab lifecycle to shutdown."""
    del request_token
    return BROWSER_LIFECYCLE_SCRIPT


DELETE_SPECIES_MODAL = r'''
<div class="species-modal-backdrop hidden" id="delete_species_modal" role="dialog" aria-modal="true" aria-labelledby="delete_species_title">
<section class="species-delete-card">
<button type="button" class="welcome-close" id="close_delete_species" aria-label="Close remove species dialog">×</button>
<h2 id="delete_species_title">Remove index or custom species?</h2>
<p id="delete_species_text">The exact cleanup list will be checked before deletion. Built-in presets and reference data are protected.</p>
<p class="delete-assets" id="delete_species_assets" aria-live="polite"></p>
<label class="delete-confirm"><input type="checkbox" id="confirm_delete_species"><span>I understand what will be removed and want to continue.</span></label>
<div class="error delete-message" id="delete_species_message"></div>
<div class="welcome-actions"><button type="button" class="secondary" id="cancel_delete_species">Cancel</button><button type="button" id="confirm_delete_species_button" disabled>Continue</button></div>
</section>
</div>'''


INDEX_ANNOTATION_MODAL = r'''
<div class="species-modal-backdrop hidden" id="index_annotation_modal" role="dialog" aria-modal="true" aria-labelledby="index_annotation_title">
<section class="welcome-card">
<button type="button" class="welcome-close" id="close_index_annotation" aria-label="Close annotation database choice">×</button>
<h2 id="index_annotation_title">Prepare this index for intron design?</h2>
<p id="index_annotation_species">The genome index will be built for the selected assembly.</p>
<p>Pre-mRNA design needs a searchable genomic annotation database that links RefSeq transcripts to their exons and introns. Building it now scans the matching GFF3 annotation and may take several minutes, but it makes future intron designs much faster.</p>
<label class="delete-confirm"><input type="checkbox" id="index_annotation_database_choice"><span>Build the genomic annotation database now for future Pre-mRNA/intron designs.</span></label>
<p class="small-note">Leave this unchecked if you will design only mature transcripts. The genome index will still be built. If you later choose Pre-mRNA, HCRProbeForge will build the annotation database during that first Pre-mRNA run.</p>
<div class="welcome-actions"><button type="button" class="secondary" id="cancel_index_annotation">Cancel</button><button type="button" id="continue_index_build">Continue index build</button></div>
</section>
</div>'''


def _render_form(request_token: str | None = None) -> str:
    request_token = request_token or secrets.token_urlsafe(32)
    supported = references.supported_species()
    # Enumerate the registry once. The old implementation independently
    # rescanned reference.json files for each preset and then rescanned them
    # again while rendering alternate species options. That was especially
    # noticeable immediately after a genome-index build, when the user was
    # sent back to this page.
    installed_references = references.list_installed_references()
    metadata_by_alias: dict[str, dict[str, object]] = {}
    for row in installed_references:
        if row.get("status") != "ready":
            continue
        for value in (
            row.get("index_species_alias"),
            row.get("species"),
            row.get("display_name"),
            row.get("scientific_name"),
        ):
            key = references.canonical_species(str(value or ""))
            if key:
                metadata_by_alias[key] = row
    # The option selected in the setup form may use a custom preset key while
    # the reference metadata uses an assembly alias. Add a direct preset-key
    # mapping so the initial badge does not depend on a later browser refresh.
    for preset in supported:
        preset_keys = {
            references.canonical_species(value)
            for value in (
                preset.key,
                preset.display_name,
                preset.scientific_name,
                *preset.aliases,
            )
            if value
        }
        matching = next(
            (
                row
                for row in installed_references
                if row.get("status") == "ready"
                and preset_keys.intersection(
                    {
                        references.canonical_species(str(value))
                        for value in (
                            row.get("index_species_alias"),
                            row.get("species"),
                            row.get("display_name"),
                            row.get("scientific_name"),
                        )
                        if value
                    }
                )
            ),
            None,
        )
        if matching is not None:
            metadata_by_alias[references.canonical_species(preset.key)] = matching
    ready = {
        preset.key: str(
            references.registered_index_is_ready(
                preset.key,
                metadata=metadata_by_alias.get(references.canonical_species(preset.key)),
            )
        ).lower()
        for preset in supported
    }
    annotation_ready = {
        preset.key: str(
            _annotation_database_ready(
                preset.key,
                metadata=metadata_by_alias.get(references.canonical_species(preset.key)),
                structural_check=False,
            )
        ).lower()
        for preset in supported
    }
    form = (
        FORM.replace("__STYLE__", STYLE)
        .replace("__FAVICON__", FAVICON_TAG)
        .replace("__VERSION__", core.__version__)
        .replace("__LOGO__", LOGO_SVG)
        .replace("__XTR_INDEX_READY__", ready["xtr"])
        .replace("__XLA_INDEX_READY__", ready["xla"])
        .replace("__ZEBRAFISH_INDEX_READY__", ready["zebrafish"])
        .replace("__MOUSE_INDEX_READY__", ready["mouse"])
        .replace("__CHICKEN_INDEX_READY__", ready["chicken"])
        .replace("__HUMAN_INDEX_READY__", ready["human"])
        .replace("__XTR_ANNOTATION_READY__", annotation_ready["xtr"])
        .replace("__XLA_ANNOTATION_READY__", annotation_ready["xla"])
        .replace("__ZEBRAFISH_ANNOTATION_READY__", annotation_ready["zebrafish"])
        .replace("__MOUSE_ANNOTATION_READY__", annotation_ready["mouse"])
        .replace("__CHICKEN_ANNOTATION_READY__", annotation_ready["chicken"])
        .replace("__HUMAN_ANNOTATION_READY__", annotation_ready["human"])
        .replace(
            "__CUSTOM_SPECIES_OPTIONS__",
            _species_options_html(
                index_ready={key: value == "true" for key, value in ready.items()},
                annotation_ready={key: value == "true" for key, value in annotation_ready.items()},
                installed_references=installed_references,
                reference_metadata=metadata_by_alias,
                annotation_structural_check=False,
            ),
        )
    )
    form = form.replace(
        "Required fields define transcript lookup and the user-facing name. Assembly details and a local genome FASTA are optional; you can supply them now or later in the Build a genome index workflow.",
        "Required fields define transcript lookup and the user-facing name. To build an index, provide either a versioned NCBI assembly accession or a local genome FASTA; these fields can be left blank only when saving the preset for later.",
    )
    form = form.replace(
        'Assembly accession <span class="optional">(optional)</span>',
        'NCBI assembly accession <span class="required">*</span> <span class="conditional">(required unless a local FASTA is supplied)</span>',
    )
    form = form.replace(
        "A stable short identifier used by the CLI, index registry, cache, and output folders.",
        "A stable technical alias used by the CLI and index registry. The readable display name is used in the webapp and reports.",
    )
    form = form.replace(
        'placeholder="e.g. versioned GCF_ or GCA_ accession"',
        'placeholder="Required for an NCBI download"',
    )
    form = form.replace(
        "builds the Bowtie2 index and registers the reference. The optional annotation database is created only when requested in the confirmation dialog; otherwise it is deferred until the first Pre-mRNA design. Uncompressed reference files are retained for reuse, and duplicate temporary <code>.gz</code> files are removed.",
        "builds the Bowtie2 index, and registers the reference. The optional annotation database is created only when requested in the confirmation dialog; otherwise it is deferred until the first Pre-mRNA design. Uncompressed reference files are retained for reuse, and duplicate temporary <code>.gz</code> files are removed.",
    )
    form = form.replace(
        'Use this when the genome is local or unavailable from NCBI. Transcript lookup still uses the NCBI organism above.',
        'Use this as the alternative genome source when you are not downloading from NCBI. Transcript lookup still uses the NCBI organism above.',
    )
    form = form.replace("option.dataset.indexReady='false';", "option.dataset.custom='true';option.dataset.localFasta=String(!!data.local_genome_fasta);option.dataset.indexReady='false';option.dataset.annotationReady='false';")
    form = form.replace("</main><script>", INDEX_ANNOTATION_MODAL + DELETE_SPECIES_MODAL + "</main><script>")
    enhancements = WEBAPP_ENHANCEMENTS.replace(
        "const removeSpeciesButton=",
        f"const hcrRequestToken={json.dumps(request_token)};const removeSpeciesButton=",
    ).replace(
        "fetch('/open-index-folder',{method:'POST'});",
        f"fetch('/open-index-folder',{{method:'POST',headers:{{'{WEB_TOKEN_HEADER}':hcrRequestToken}}}});",
    ).replace(
        "const payload=new Blob([''],{type:'text/plain'});",
        "const payload=new Blob([hcrRequestToken],{type:'text/plain'});",
    )
    # Built-in presets remain in the selector after their installed index is
    # removed; custom presets retain the existing full-removal behavior.
    enhancements = enhancements.replace(
        "removeSpeciesButton.hidden=!custom;",
        "removeSpeciesButton.hidden=!(custom||ready);",
    ).replace(
        "removeSpeciesButton.dataset.help=custom?'Remove this custom species preset, its saved HCRProbeForge data, and its unshared Bowtie2 index files. Built-in presets are protected.':'Remove this custom species preset.';",
        "removeSpeciesButton.dataset.help=custom?'Remove this custom species preset, its saved HCRProbeForge data, and its unshared Bowtie2 index files.':ready?'Remove only this built-in species index files and registration. The built-in preset, genome, and annotation data are retained.':'Remove the installed index for this species.';",
    ).replace(
        "if(option&&option.dataset.custom==='true')openDeleteSpecies(option);",
        "if(option&&(option.dataset.custom==='true'||option.dataset.indexReady==='true'))openDeleteSpecies(option);",
    ).replace(
        "option.remove();species.selectedIndex=0;",
        "if(data.builtin!=='true'){option.remove();species.selectedIndex=0;}else{option.dataset.indexReady='false';}",
    )
    enhancements = enhancements.replace(
        "deleteSpeciesText.textContent=`This removes “${option.textContent}” from the species list, its saved preset data, stale config registrations, and the listed unshared index files. Built-in presets are protected.`;",
        "deleteSpeciesText.textContent=plan.builtin?`This keeps the built-in “${option.textContent}” preset and all genome/annotation reference files. It removes only the listed unshared Bowtie2 index files and matching index registrations.`:`This removes “${option.textContent}” from the species list, its saved preset data, stale config registrations, and the listed unshared index files.`;",
    )
    return form.replace("</script></body></html>", enhancements + "</script></body></html>")


def _folder_component(value: str) -> str:
    return core.safe_name(value.strip()) or "workflow"


def _field_enabled(fields: dict[str, str], name: str) -> bool:
    """Interpret checkbox values from both HTML forms and CLI-like callers."""
    return str(fields.get(name, "")).strip().casefold() in {"1", "true", "yes", "on"}


def _workflow_root(fields: dict[str, str], mode: str, base: Path, input_path: Path | None = None) -> Path:
    """Return a persistent, species/workflow-specific output directory."""
    species = fields.get("species", "xtr").strip() or "xtr"
    requested = fields.get("folder_name", "").strip()
    if mode == "manifest":
        return references.results_root(base)
    workflow_root = references.species_run_root(base, species, mode)
    if requested:
        return workflow_root / _folder_component(requested)
    if mode == "plot" and input_path is not None:
        return workflow_root / _folder_component(f"{input_path.stem}_plot")
    if mode == "qc" and input_path is not None:
        return workflow_root / _folder_component(f"{input_path.stem}_qc")
    return workflow_root


_RUN_OUTPUT_SUFFIXES = (
    "_final_selected_pairs.tsv",
    "_final_IDT_order.csv",
    "_final_probe_map.png",
    "_final_probe_map.svg",
    "_final_oligo_structure_QC.xlsx",
    "_selected_pairs_probe_map.png",
    "_selected_pairs_probe_map.svg",
    "_probe_map.png",
    "_probe_map.svg",
    "_oligo_structure_QC.xlsx",
    "_auto_curate_report.json",
    "_auto_curate_report.md",
    "_summary.json",
)


def _result_target_roots(console_text: str) -> list[Path]:
    """Extract target directories named by this run's own output messages."""
    roots: list[Path] = []
    seen: set[Path] = set()
    for line in console_text.splitlines():
        if ":" not in line:
            continue
        # The first colon separates the console label from the path.  The
        # remainder may itself contain a Windows drive colon (``C:\\...``),
        # so never split the complete line on every colon.
        candidate = Path(line.partition(":")[2].strip().strip('"')).expanduser()
        if not candidate.name.endswith(_RUN_OUTPUT_SUFFIXES):
            continue
        root = candidate.parent
        if root.name == "details":
            root = root.parent
        try:
            root = root.resolve()
        except OSError:
            continue
        if root.is_dir() and root not in seen:
            roots.append(root)
            seen.add(root)
    return roots


def _common_directory(paths: list[Path], fallback: Path) -> Path:
    if not paths:
        return fallback.resolve()
    try:
        return Path(os.path.commonpath([str(path.resolve()) for path in paths]))
    except (OSError, ValueError):
        return fallback.resolve()


def _root_context(root: Path, output_directory: str) -> str:
    """Return the useful target/channel context for a result root."""
    try:
        relative = root.resolve().relative_to(Path(output_directory).expanduser().resolve())
    except (OSError, ValueError):
        relative = Path(root.name)
    value = relative.as_posix()
    return value if value not in {"", "."} else root.name


def _display_file_path(path: Path, root: Path, output_directory: str) -> str:
    context = _root_context(root, output_directory)
    return f"{context}/{path.name}" if context else path.name


def _map_display_label(path: Path, root: Path, output_directory: str) -> str:
    """Use target and channel context so repeated transcript filenames differ."""
    context = _root_context(root, output_directory)
    components = list(root.parts)
    channel = next(
        (component for component in reversed(components) if re.search(r"(?:^|_)B[1-5]$", component)),
        "",
    )
    target = root.name
    if channel and channel not in context:
        return f"{channel} · {target} · {path.name}"
    return f"{context} · {path.name}" if context else path.name


def _web_result(
    code: int,
    stdout: str,
    stderr: str,
    roots: list[Path],
    **extra: object,
) -> dict[str, object]:
    unique_roots: list[str] = []
    for root in roots:
        resolved = str(Path(root).expanduser().resolve())
        if resolved not in unique_roots:
            unique_roots.append(resolved)
    return {
        "return_code": code,
        "success": code == 0,
        "stdout": stdout,
        "stderr": stderr,
        "roots": unique_roots,
        "output_directory": unique_roots[0] if unique_roots else "",
        **extra,
    }


def _run_submission(
    fields: dict[str, str],
    uploads: list[tuple[str, str, bytes]],
    progress_callback: ProgressCallback | None = None,
    cancel_event: threading.Event | None = None,
) -> dict[str, object]:
    check = cancel_event.is_set if cancel_event is not None else None
    with core.cancellation_scope(check):
        core.raise_if_cancelled()
        return _run_submission_inner(fields, uploads, progress_callback=progress_callback)


def _run_submission_inner(
    fields: dict[str, str],
    uploads: list[tuple[str, str, bytes]],
    progress_callback: ProgressCallback | None = None,
) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="hcrprobeforge-web-") as temporary:
        return _run_submission_in_directory(
            fields,
            uploads,
            progress_callback=progress_callback,
            job_dir=Path(temporary),
        )


def _run_submission_in_directory(
    fields: dict[str, str],
    uploads: list[tuple[str, str, bytes]],
    progress_callback: ProgressCallback | None = None,
    *,
    job_dir: Path,
) -> dict[str, object]:
    mode = fields.get("mode", "design").strip().lower()
    _validate_form_values(fields)
    if fields.get("transcript_policy", "auto").strip() == "interactive":
        raise ValueError("interactive transcript selection is available from the CLI; choose auto, longest, or require-accession in the webapp")
    output_base = _output_path(fields)
    if mode == "index":
        species = fields.get("species", "").strip()
        preset = references.get_species_preset(species)
        if preset is None:
            raise ValueError("Choose a built-in or user-created species preset before building an index")
        try:
            threads = int(fields.get("index_threads", "4"))
        except ValueError as exc:
            raise ValueError("Build threads must be a whole number") from exc
        if threads < 1:
            raise ValueError("Build threads must be at least 1")
        force = _field_enabled(fields, "index_force")
        build_annotation_database = _field_enabled(fields, "index_annotation_database")
        updates: list[str] = []

        def index_progress(update: dict[str, object]) -> None:
            if progress_callback is not None:
                progress_callback(update)

        local_fasta = _upload_or_path(fields, uploads, "index_fasta_path", "index_fasta_upload", job_dir)
        local_annotation = _upload_or_path(
            fields,
            uploads,
            "index_annotation_path",
            "index_annotation_upload",
            job_dir,
        )
        local_fasta = _preserve_uploaded_input(
            local_fasta,
            uploads,
            "index_fasta_upload",
            output_base,
        )
        if local_fasta is None and preset.local_genome_fasta:
            saved_fasta = Path(preset.local_genome_fasta).expanduser()
            if saved_fasta.is_file():
                local_fasta = saved_fasta
        assembly_accession = fields.get("index_assembly_accession", "").strip() or None
        assembly_name = fields.get("index_assembly_name", "").strip() or None
        index_alias = fields.get("index_alias", "").strip() or None
        if local_fasta is not None:
            result = references.build_local_index(
                preset.key,
                local_fasta,
                threads=threads,
                force=force,
                assembly=assembly_name or preset.assembly_name,
                annotation=local_annotation,
                build_annotation_database=build_annotation_database,
                progress_callback=index_progress,
            )
        else:
            if not assembly_accession and not preset.assembly_accession:
                raise ValueError(
                    f"Custom preset {preset.display_name} was saved, but no genome source is configured yet. "
                    "In Build a genome index, provide a versioned NCBI assembly accession or upload a local genome FASTA."
                )
            result = references.fetch_and_build_index(
                preset.key,
                assembly_accession=assembly_accession,
                assembly_name=assembly_name,
                index_alias=index_alias,
                threads=threads,
                force=force,
                build_annotation_database=build_annotation_database,
                email=fields.get("email") or os.getenv("NCBI_EMAIL"),
                api_key=fields.get("api_key") or os.getenv("NCBI_API_KEY"),
                progress_callback=index_progress,
            )
        metadata_path = Path(str(result.get("metadata_path") or references.species_data_root()))
        alias = core.safe_name(str(result.get("index_species_alias") or preset.key)) or preset.key
        # The HCRProbeDesign index itself stays in its persistent data
        # directory.  Only a clearly labelled project record is written next
        # to hcr_results, so result folders never look like genome indexes.
        assembly_component = core.safe_name(str(result.get("assembly") or result.get("assembly_accession") or "local_reference")) or "local_reference"
        index_output = references.index_metadata_root(output_base) / references.species_component(preset.key) / assembly_component
        index_output.mkdir(parents=True, exist_ok=True)
        web_metadata = index_output / "reference.json"
        index_prefix = Path(str(result.get("index_prefix") or references.registered_index_prefix(alias))).expanduser().resolve()
        result["index_prefix"] = str(index_prefix)
        result["web_metadata_path"] = str(web_metadata)
        web_metadata.write_text(json.dumps(core.json_ready(result), indent=2, sort_keys=True) + "\n", encoding="utf-8")
        stdout = (
            f"Reference ready: {preset.display_name} · {result.get('assembly')} · "
            f"{result.get('assembly_accession')}\n"
            f"Bowtie2 index prefix: {index_prefix}\n"
            f"Index metadata: {web_metadata}\n"
            f"HCRProbeDesign registry metadata: {metadata_path}\n"
            "The genome was passed to buildGenomeIndex and registered for future designs.\n"
        )
        return _web_result(
            0,
            stdout,
            "",
            [index_output],
            workflow="index",
            output_directory=str(index_output),
            reference=result,
            index_prefix=str(index_prefix),
            index_metadata_path=str(web_metadata),
        )
    if mode == "design":
        gene = fields.get("gene", "").strip() or None
        accession = fields.get("accession", "").strip() or None
        fasta = _upload_or_path(fields, uploads, "fasta_path", "fasta_upload", job_dir)
        if fasta and (gene or accession):
            raise ValueError("FASTA input cannot be combined with a gene symbol or accession")
        if not fasta and not (gene or accession):
            raise ValueError("Provide a gene, accession, or FASTA input")
        output = _workflow_root(fields, mode, output_base)
        fasta = _preserve_uploaded_input(fasta, uploads, "fasta_upload", output)
        forwarded = _common_design_args(fields)
        if "--_species-scoped-outdir" not in forwarded:
            forwarded.append("--_species-scoped-outdir")
        channel = fields.get("channel", "ALL").strip().upper()
        if channel == "ALL":
            result = batch.run_all_channels(
                gene=gene,
                accession=accession,
                fasta=fasta,
                output_root=output,
                extra_args=forwarded,
                capture=True,
                smart_default=False,
                progress_callback=progress_callback,
                output_scoped=True,
                workflow="design",
                cache_root=references.species_cache_root(output_base, fields.get("species", "xtr"), "design"),
            )
            console = str(result["stdout"]) + str(result["stderr"])
            roots = _result_target_roots(console)
            if not roots and result.get("success"):
                roots = [Path(str(record["output_dir"])) for record in result.get("channels", [])]
            reported_roots = roots if roots else (
                [output] if int(result.get("return_code", 1)) == 0 else []
            )
            return _web_result(
                int(result["return_code"]),
                str(result["stdout"]),
                str(result["stderr"]),
                reported_roots,
                workflow="design",
                output_directory=(
                    str(_common_directory(reported_roots, output / _folder_component(gene or accession or (fasta.stem if fasta else "custom_target"))))
                    if reported_roots else ""
                ),
            )
        argv: list[str] = []
        if gene:
            argv.append(gene)
        if accession:
            argv += ["--accession", accession]
        if fasta:
            argv += ["--fasta", str(fasta)]
        argv += ["--channel", channel, "--outdir", str(output), "--_species-scoped-outdir"]
        argv += forwarded
        # This is the first notification for a single-target run. It must
        # enter the normal starting range; announcing the later design phase
        # here made a one-channel run appear to begin at 60%.
        _notify_web_progress(progress_callback, phase="starting", message="Starting HCRProbeForge design")
        code, stdout, stderr = _invoke_with_cache(
            argv,
            output_base=output_base,
            species=fields.get("species", "xtr"),
            workflow="design",
            progress_callback=progress_callback,
        )
        roots = _result_target_roots(stdout + stderr)
        reported_roots = roots if roots else ([output] if code == 0 else [])
        return _web_result(
            code,
            stdout,
            stderr,
            reported_roots,
            workflow="design",
            output_directory=str(_common_directory(reported_roots, output)) if reported_roots else "",
        )
    if mode == "manifest":
        manifest = _upload_or_path(fields, uploads, "manifest_path", "manifest_upload", job_dir)
        if not manifest:
            raise ValueError("Upload a manifest or provide its local path")
        output = _workflow_root(fields, mode, output_base)
        manifest = _preserve_uploaded_input(manifest, uploads, "manifest_upload", output)
        output.mkdir(parents=True, exist_ok=True)
        channel = fields.get("channel", "ALL").strip().upper()
        channels = batch.CHANNELS if channel == "ALL" else (channel,)
        result = batch.run_manifest(
            manifest=manifest,
            project_root=output,
            extra_args=_common_design_args(fields),
            channels=channels,
            smart_default=False,
            progress_callback=progress_callback,
            run_name=fields.get("folder_name", "").strip() or manifest.stem,
        )
        manifest_root = Path(str(result.get("project_root") or output / "manifest"))
        roots = [Path(str(root)).expanduser().resolve() for root in result.get("roots", [])]
        if not roots:
            roots = _result_target_roots(str(result.get("stdout", "")) + str(result.get("stderr", "")))
        return _web_result(
            int(result["return_code"]),
            str(result["stdout"]),
            str(result["stderr"]),
            roots,
            workflow="manifest",
            output_directory=str(manifest_root),
        )
    existing = _upload_or_path(fields, uploads, "existing_path", "existing_upload", job_dir)
    if not existing:
        raise ValueError("Upload an input file or provide its local path")
    output = _workflow_root(fields, mode, output_base, existing)
    existing = _preserve_uploaded_input(existing, uploads, "existing_upload", output)
    if mode == "plot":
        argv = ["--plot-only", str(existing), "--outdir", str(output), "--species", fields.get("species", "xtr"), "--_species-scoped-outdir"] + _plot_args(fields)
        _notify_web_progress(progress_callback, phase="plot", message="Creating the probe map")
        code, stdout, stderr = _invoke_with_cache(
            argv,
            output_base=output_base,
            species=fields.get("species", "xtr"),
            workflow="plot",
        )
        return _web_result(code, stdout, stderr, [output] if code == 0 else [], workflow="plot")
    if mode == "qc":
        argv = ["--qc-only", str(existing)]
        qc_value = fields.get("qc_output", "").strip()
        qc_output = Path(qc_value).expanduser() if qc_value else output / f"{_folder_component(existing.stem)}_oligo_structure_QC.xlsx"
        if not qc_output.is_absolute():
            qc_output = output / qc_output
        qc_output = qc_output.resolve()
        argv += ["--qc-output", str(qc_output)]
        _notify_web_progress(progress_callback, phase="qc", message="Running oligo QC and writing the workbook")
        code, stdout, stderr = _invoke_with_cache(
            argv,
            output_base=output_base,
            species=fields.get("species", "xtr"),
            workflow="qc",
        )
        return _web_result(code, stdout, stderr, [qc_output.parent] if code == 0 else [], workflow="qc")
    raise ValueError(f"unknown workflow: {mode}")


def _notify_web_progress(progress_callback: ProgressCallback | None, **update: object) -> None:
    if progress_callback is None:
        return
    try:
        progress_callback(update)
    except core.RunCancelled:
        raise
    except Exception:
        pass


def _validate_form_values(fields: dict[str, str]) -> None:
    """Validate browser values before a worker is started."""
    mode = fields.get("mode", "design").strip().lower()
    if mode not in {"design", "manifest", "plot", "qc", "index"}:
        raise ValueError("Choose one of the available workflows")
    channel = fields.get("channel", "ALL").strip().upper()
    if channel != "ALL" and channel not in batch.CHANNELS:
        raise ValueError("Channel must be B1, B2, B3, B4, B5, or All channels")
    integer_fields = {
        "tile_size": 1,
        "max_run_mismatches": 0,
        "max_probes": 1,
        "num_hits_allowed": 1,
        "plot_dpi": 72,
        "transcript_length": 1,
        "target_probes": 1,
        "min_acceptable_probes": 1,
        "max_auto_runs": 0,
        "auto_curate_max_probes": 1,
        "index_threads": 1,
    }
    for field, minimum in integer_fields.items():
        value = fields.get(field, "").strip()
        if not value:
            continue
        try:
            parsed = int(value)
        except ValueError as exc:
            raise ValueError(f"{field.replace('_', ' ').capitalize()} must be a whole number") from exc
        if parsed < minimum:
            raise ValueError(f"{field.replace('_', ' ').capitalize()} must be at least {minimum}")
    float_fields = ("min_gc", "max_gc", "min_gibbs", "max_gibbs", "target_gibbs", "dtm_max")
    parsed_float: dict[str, float] = {}
    for field in float_fields:
        value = fields.get(field, "").strip()
        if not value:
            continue
        try:
            parsed_float[field] = float(value)
        except ValueError as exc:
            raise ValueError(f"{field.replace('_', ' ').capitalize()} must be a number") from exc
    if "min_gc" in parsed_float and "max_gc" in parsed_float and not 0 <= parsed_float["min_gc"] < parsed_float["max_gc"] <= 100:
        raise ValueError("Minimum GC must be lower than maximum GC and both must be between 0 and 100")
    if "min_gibbs" in parsed_float and "max_gibbs" in parsed_float and parsed_float["min_gibbs"] >= parsed_float["max_gibbs"]:
        raise ValueError("Minimum Gibbs must be lower than maximum Gibbs")
    if "target_gibbs" in parsed_float and {"min_gibbs", "max_gibbs"}.issubset(parsed_float):
        if not parsed_float["min_gibbs"] <= parsed_float["target_gibbs"] <= parsed_float["max_gibbs"]:
            raise ValueError("Target Gibbs must lie between the minimum and maximum Gibbs values")
    if "dtm_max" in parsed_float and parsed_float["dtm_max"] < 0:
        raise ValueError("Maximum dTm cannot be negative")
    for field, allowed in {
        "plot_theme": {"pastel", "minimal"},
        "plot_color_by": {"gc", "dtm", "order"},
        "transcript_policy": {"auto", "longest", "interactive", "require-accession"},
        "target_type": {"mature", "pre-mrna", "pre-mrna-whole"},
        "qc_stringency": {"strict", "balanced", "permissive"},
        "auto_curate_plan": {"standard", "conservative", "deep"},
        "coverage_policy": {"balanced", "qc-only"},
    }.items():
        value = fields.get(field, "").strip()
        if value and value not in allowed:
            raise ValueError(f"Unsupported {field.replace('_', ' ')} value: {value}")
    target_type = fields.get("target_type", "mature").strip() or "mature"
    if fields.get("premrna_introns", "").strip():
        if target_type != "pre-mrna":
            raise ValueError("Intron numbers require the Pre-mRNA target type")
        from .premrna import parse_intron_selection

        parse_intron_selection(fields.get("premrna_introns"))
    if target_type in {"pre-mrna", "pre-mrna-whole"} and mode != "design":
        raise ValueError("Pre-mRNA target selection is currently available for Design one target only")
    if mode == "index":
        accession = fields.get("index_assembly_accession", "").strip()
        if accession and not references.ASSEMBLY_ACCESSION_RE.fullmatch(accession):
            raise ValueError("Assembly accession must look like GCF_<ASSEMBLY_ID>.<VERSION> or GCA_<ASSEMBLY_ID>.<VERSION>")
        alias = fields.get("index_alias", "").strip()
        if alias and core.safe_name(alias) != alias:
            raise ValueError("Index alias may contain only letters, numbers, dots, underscores, and hyphens")
        for field in ("index_assembly_name", "index_alias"):
            if any(character in fields.get(field, "") for character in "\r\n\t"):
                raise ValueError(f"{field.replace('_', ' ').capitalize()} cannot contain line breaks or tabs")
    qc_output = fields.get("qc_output", "").strip()
    if qc_output and Path(qc_output).suffix.lower() not in {".xlsx", ".xlsm"}:
        raise ValueError("QC workbook output must end in .xlsx or .xlsm")
    if "target_probes" in fields and "min_acceptable_probes" in fields:
        target_value = fields.get("target_probes", "").strip()
        minimum_value = fields.get("min_acceptable_probes", "").strip()
        if target_value and minimum_value and int(minimum_value) > int(target_value):
            raise ValueError("Minimum acceptable pairs cannot be greater than the target probe pairs")

    # The browser always submits the species select.  Keep the historical
    # programmatic API forgiving for plot/QC callers that omit it: those
    # workflows do not use a species to resolve an NCBI transcript.
    species = fields.get("species", "xtr").strip() or "xtr"
    if not species:
        raise ValueError("Choose a species before continuing")
    preset = references.get_species_preset(species)
    installed = references.find_installed_reference(species) if preset is None else None
    if mode == "index" and preset is None:
        raise ValueError("Choose a built-in or user-created species preset before building an index")
    if preset is None and installed is None:
        raise ValueError("Choose a supported species or a registered genome assembly")
    expected_organism = (
        preset.scientific_name
        if preset is not None
        else str(installed.get("scientific_name") or "") if installed is not None else ""
    ).strip()
    supplied_organism = fields.get("organism", "").strip()
    if expected_organism and supplied_organism and supplied_organism.casefold() != expected_organism.casefold():
        raise ValueError(
            f"NCBI organism must match the selected species ({expected_organism}); "
            f"received {supplied_organism}."
        )
    if mode in {"design", "manifest"} and "no_genomemask" not in fields and not fields.get("index", "").strip():
        if not references.registered_index_is_ready(species):
            raise ValueError(core.missing_genome_index_message(species))
    if mode in {"design", "manifest"}:
        if "auto_curate" in fields and "no_oligo_qc" in fields:
            raise ValueError("Smart auto-curation cannot be combined with Skip final oligo QC")


def _inspect_submission(fields: dict[str, str], uploads: list[tuple[str, str, bytes]]) -> dict[str, str]:
    with tempfile.TemporaryDirectory(prefix="hcrprobeforge-inspect-") as temporary:
        return _inspect_submission_in_directory(fields, uploads, Path(temporary))


def _inspect_submission_in_directory(
    fields: dict[str, str],
    uploads: list[tuple[str, str, bytes]],
    job_dir: Path,
) -> dict[str, str]:
    mode = fields.get("mode", "design").strip().lower()
    _validate_form_values(fields)
    transcript_policy = fields.get("transcript_policy", "auto").strip()
    if mode in {"design", "manifest"} and transcript_policy == "interactive":
        raise ValueError("Interactive transcript selection is available in the CLI; choose auto, longest, or require-accession here")
    if mode == "index":
        preset = references.get_species_preset(fields.get("species", ""))
        if preset is None:
            raise ValueError("Choose a built-in or user-created species preset before building an index")
        try:
            threads = int(fields.get("index_threads", "4"))
        except ValueError as exc:
            raise ValueError("Build threads must be a whole number") from exc
        if threads < 1:
            raise ValueError("Build threads must be at least 1")
        local_fasta = _upload_or_path(fields, uploads, "index_fasta_path", "index_fasta_upload", job_dir)
        local_annotation = _upload_or_path(
            fields,
            uploads,
            "index_annotation_path",
            "index_annotation_upload",
            job_dir,
        )
        if local_fasta is None and preset.local_genome_fasta:
            saved_fasta = Path(preset.local_genome_fasta).expanduser()
            if saved_fasta.is_file():
                local_fasta = saved_fasta
        accession = fields.get("index_assembly_accession", "").strip() or preset.assembly_accession
        assembly = fields.get("index_assembly_name", "").strip() or preset.assembly_name
        alias = fields.get("index_alias", "").strip() or "automatic alias for this assembly"
        annotation_database_requested = _field_enabled(fields, "index_annotation_database")
        if local_fasta is not None:
            if annotation_database_requested and local_annotation is None:
                raise ValueError(
                    "The annotation database option is enabled, but no matching local GFF3 was supplied. "
                    "Provide the GFF3 from the same assembly or leave the option unchecked."
                )
            annotation_note = (
                f" with matching genomic annotation {local_annotation.name}"
                if local_annotation is not None
                else " without a local annotation"
            )
            database_note = (
                " The optional annotation database will also be built."
                if annotation_database_requested
                else " The annotation database will be deferred until the first Pre-mRNA run."
            )
            return {
                "message": (
                    f"Ready to build {preset.display_name} from local genome FASTA {local_fasta.name} "
                    f"as {fields.get('index_alias', '').strip() or preset.key}{annotation_note}, "
                    f"using {threads} thread(s).{database_note}"
                )
            }
        if not accession:
            raise ValueError(
                f"Custom preset {preset.display_name} has no genome source yet. "
                "Provide a versioned NCBI assembly accession or upload a local genome FASTA before building the index."
            )
        if not references.ASSEMBLY_ACCESSION_RE.fullmatch(accession):
            raise ValueError("Assembly accession must look like GCF_<ASSEMBLY_ID>.<VERSION> or GCA_<ASSEMBLY_ID>.<VERSION>")
        return {
            "message": (
                f"Ready to build {preset.display_name}: {assembly} ({accession}) as {alias}, "
                f"using {threads} thread(s). NCBI will be checked for the exact assembly "
                "before the genome is downloaded and indexed. The annotation database will "
                f"{'also be built now' if annotation_database_requested else 'be deferred until the first Pre-mRNA run'}."
            )
        }
    if mode == "design":
        design_accession = fields.get("accession", "").strip()
        if transcript_policy == "require-accession" and not design_accession:
            raise ValueError("Transcript policy 'require-accession' needs a RefSeq transcript accession")
    if mode == "manifest":
        path = _upload_or_path(fields, uploads, "manifest_path", "manifest_upload", job_dir)
        if not path:
            raise ValueError("Upload a manifest or provide its local path")
        if path.suffix.lower() not in {".tsv", ".csv", ".xlsx", ".xlsm", ".txt"}:
            raise ValueError("Manifest must be CSV, TSV, XLSX, XLSM, or text")
        all_rows, collisions = batch._preflight_all_rows(path)
        jobs = [row for row in all_rows if row.get("status") in {"valid", "shared_accession"}]
        if not all_rows:
            raise ValueError("Manifest is empty. Add at least one gene symbol or accession row.")
        if collisions:
            details = "; ".join(
                f"line {row.get('line_number')}: {row.get('gene_symbol')} -> {row.get('gene_dir')}"
                for row in collisions
            )
            raise ValueError(f"Manifest has directory-name collisions: {details}")
        global_channel = fields.get("channel", "ALL").strip().upper()
        manifest_channels = batch.CHANNELS if global_channel == "ALL" else (global_channel,)
        channel_collisions = batch.manifest_channel_collisions(all_rows, manifest_channels)
        if channel_collisions:
            details = "; ".join(
                f"line {row.get('line_number')}: {row.get('gene_symbol')} ({row.get('channel') or 'global channels'})"
                for row in channel_collisions
            )
            raise ValueError(
                "Manifest assigns the same target/channel output more than once: "
                f"{details}. Use one row per target/channel."
            )
        if not jobs:
            raise ValueError(
                "Manifest has no structurally valid rows. Check the gene_symbol, accession, "
                "and channel columns before running it."
            )
        if transcript_policy == "require-accession":
            missing_accession = sum(not str(row.get("accession") or "").strip() for row in jobs)
            if missing_accession:
                raise ValueError(
                    f"Transcript policy 'require-accession' needs an accession for every valid manifest row; "
                    f"{missing_accession} row(s) are missing one."
                )
        preflight_skips = sum(row.get("status") not in {"valid", "shared_accession"} for row in all_rows)
        per_row_channels = sum(bool(str(row.get("channel") or "").strip()) for row in all_rows)
        normalized = sum(
            str(row.get("gene_symbol")) != core.safe_name(str(row.get("gene_symbol", "")))
            for row in batch._manifest_rows(path)
        )
        channel_note = f" {per_row_channels} row(s) specify their own HCR channel." if per_row_channels else ""
        skip_note = f" {preflight_skips} structurally invalid row(s) will not run." if preflight_skips else ""
        return {
            "message": (
                f"Manifest ready: {len(jobs)} valid unique job(s);{channel_note}{skip_note} "
                f"{normalized} normalized gene label(s). NCBI will check each gene or accession "
                "against the selected species during the run; true species mismatches are "
                "skipped, unresolved inputs fail, and other rows continue."
            )
        }
    if mode in {"plot", "qc"}:
        path = _upload_or_path(fields, uploads, "existing_path", "existing_upload", job_dir)
        if not path:
            raise ValueError("Upload an input file or provide its local path")
        if mode == "plot":
            if path.suffix.lower() not in {".tsv", ".txt"}:
                raise ValueError("Plot workflow expects an HCRProbeDesign TSV or text table")
            rows = core.read_probe_rows(path)
            if not rows:
                raise ValueError("Probe-table inspection found no usable probe rows")
            final_probe_end = max(int(row["end"]) for row in rows)
            transcript_length = fields.get("transcript_length", "").strip()
            if transcript_length and int(transcript_length) < final_probe_end:
                raise ValueError(
                    f"Transcript length ({transcript_length}) is shorter than the final probe coordinate "
                    f"({final_probe_end})."
                )
            length_note = f" Transcript length: {transcript_length} nt." if transcript_length else f" Last probe coordinate: {final_probe_end}."
            return {"message": f"Probe table ready: {path.name}; {len(rows)} probe row(s) with valid coordinates.{length_note}"}
        if path.suffix.lower() not in {".csv", ".tsv", ".txt", ".xlsx", ".xlsm"}:
            raise ValueError("QC workflow expects CSV, TSV, TXT, XLSX, or XLSM")
        table = core.read_oligo_table(path)
        if getattr(table, "empty", False):
            raise ValueError("QC inspection found no oligo rows")
        return {"message": f"Oligo table ready: {path.name}; {len(table)} oligo row(s) detected."}
    path = _upload_or_path(fields, uploads, "fasta_path", "fasta_upload", job_dir)
    gene = fields.get("gene", "").strip()
    accession = fields.get("accession", "").strip()
    if path and (gene or accession):
        raise ValueError("FASTA input cannot be combined with a gene symbol or accession")
    if path and "all_transcripts" in fields:
        raise ValueError("FASTA input cannot be combined with Design all linked transcripts")
    if path and transcript_policy != "auto":
        raise ValueError("Transcript selection policy applies to NCBI gene lookup, not FASTA input")
    if path:
        records = core.parse_fasta_records(path)
        if not records:
            raise ValueError("No FASTA record was found")
        if len(records) != 1:
            raise ValueError("Design input must contain exactly one FASTA record")
        lengths = ", ".join(str(len(record.get("sequence", ""))) for record in records)
        return {
            "message": (
                f"FASTA ready for {fields.get('species', '').strip()}: {len(records)} record(s), "
                f"sequence length(s): {lengths} nt. The sequence will be used as supplied."
            )
        }
    if accession and not batch.ACCESSION_RE.fullmatch(accession):
        raise ValueError("Transcript accession should look like NM_001234.1, XM_001234.2, NR_..., or XR_...")
    if "all_transcripts" in fields and accession:
        raise ValueError("Design all linked transcripts cannot be combined with an explicit accession")
    if not references.get_species_preset(fields.get("species", "")) and not references.find_installed_reference(fields.get("species", "")):
        raise ValueError("Choose a supported species or a registered genome assembly before transcript lookup")
    if gene or accession:
        return {"message": f"NCBI input ready: gene={gene or '(none)'}, accession={accession or '(automatic selection)'}."}
    raise ValueError("Provide a gene, accession, or FASTA input")


def _choose_local_folder() -> str:
    if sys.platform == "darwin":
        result = subprocess.run(["osascript", "-e", 'POSIX path of (choose folder with prompt "Choose a folder")'], capture_output=True, text=True, timeout=120, check=False)
        return result.stdout.strip() if result.returncode == 0 else ""
    try:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        selected = filedialog.askdirectory(title="Choose a folder")
        root.destroy()
        return selected
    except Exception:
        return ""


def _is_wsl() -> bool:
    """Return whether the server is running inside Windows Subsystem for Linux."""
    if not sys.platform.startswith("linux"):
        return False
    if os.environ.get("WSL_INTEROP") or os.environ.get("WSL_DISTRO_NAME"):
        return True
    try:
        release = platform.release().casefold()
    except OSError:
        release = ""
    return "microsoft" in release or "wsl" in release


def _open_browser_url(url: str) -> bool:
    """Open a local URL using the host browser, including from WSL."""
    if _is_wsl():
        explorer = shutil.which("explorer.exe")
        if explorer:
            try:
                subprocess.Popen(
                    [explorer, url],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    start_new_session=True,
                )
                return True
            except OSError:
                pass
    try:
        return bool(webbrowser.open(url))
    except Exception:
        return False


def _wsl_windows_path(path: Path) -> str:
    """Convert a WSL path to a Windows path for Explorer."""
    wslpath = shutil.which("wslpath")
    if not wslpath:
        raise RuntimeError("WSL could not find wslpath; open the folder manually from Ubuntu.")
    result = subprocess.run(
        [wslpath, "-w", str(path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode or not result.stdout.strip():
        detail = result.stderr.strip() or "path conversion failed"
        raise RuntimeError(f"WSL could not convert the output path for Windows Explorer: {detail}")
    return result.stdout.strip()


def _open_local_directory(path: Path) -> None:
    """Open an output directory in the operating system's file browser."""
    directory = path.expanduser().resolve()
    if not directory.is_dir():
        raise FileNotFoundError(f"Output directory does not exist: {directory}")
    if sys.platform == "darwin":
        subprocess.Popen(["open", str(directory)])
    elif _is_wsl():
        explorer = shutil.which("explorer.exe")
        if not explorer:
            raise RuntimeError(
                "Windows Explorer (explorer.exe) was not found. Open the folder from Windows using the WSL path manually."
            )
        subprocess.Popen(
            [explorer, _wsl_windows_path(directory)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True,
        )
    elif os.name == "nt":
        os.startfile(str(directory))  # type: ignore[attr-defined]
    else:
        opener = shutil.which("xdg-open") or shutil.which("gio")
        if not opener:
            raise RuntimeError("No desktop directory opener was found (expected xdg-open or gio)")
        command = [opener, "open", str(directory)] if Path(opener).name == "gio" else [opener, str(directory)]
        subprocess.Popen(command)


def _roots_from_job(job: dict[str, object]) -> list[Path]:
    return [Path(str(root)).expanduser().resolve() for root in job.get("roots", [])]


def _file_url(token: str, root_index: int, root: Path, path: Path) -> str:
    relative = path.relative_to(root).as_posix()
    return f"/file/{quote(token)}/{root_index}/{quote(relative)}"


def _important_files(roots: list[Path]) -> list[tuple[int, Path, str]]:
    """Return only user-facing result files, excluding technical details."""
    patterns = (
        ("reference.json", "Reference metadata"),
        ("final_selected_pairs.tsv", "Selected probe pairs"),
        ("final_IDT_order.csv", "IDT order file"),
        ("final_probe_map.png", "Final probe map"),
        ("final_probe_map.svg", "Final probe map (SVG)"),
        ("final_oligo_structure_QC.xlsx", "Final oligo QC workbook"),
        ("_probe_map.png", "Probe map"),
        ("_probe_map.svg", "Probe map (SVG)"),
        ("_oligo_structure_QC.xlsx", "Oligo QC workbook"),
    )
    matches: list[tuple[int, Path, str]] = []
    seen: set[Path] = set()
    for root_index, root in enumerate(roots):
        if not root.exists():
            continue
        for path in sorted(item for item in root.rglob("*") if item.is_file()):
            if "details" in path.relative_to(root).parts:
                continue
            for suffix, label in patterns:
                if path.name.endswith(suffix):
                    if path not in seen:
                        matches.append((root_index, path, label))
                        seen.add(path)
                    break
    return matches


def _map_files(roots: list[Path]) -> list[tuple[int, Path]]:
    return [(index, path) for index, path, label in _important_files(roots) if "map" in label.lower()]


def _selected_pair_files(roots: list[Path]) -> list[tuple[int, Path, str]]:
    return [item for item in _important_files(roots) if item[1].name.endswith("_final_selected_pairs.tsv")]


def _selected_pair_label(root: Path, output_directory: str) -> str:
    context = _root_context(root, output_directory)
    channel = ""
    for component in context.split("/"):
        match = re.search(r"(?:^|_)B([1-5])$", component)
        if match:
            channel = f"B{match.group(1)}"
            break
    parts = [part for part in (channel, root.name) if part]
    return " · ".join(parts + ["Selected probe pairs"])


def _important_file_links(token: str, roots: list[Path], output_directory: str = "", workflow: str = "design") -> str:
    items: list[str] = []
    matches = _selected_pair_files(roots) if workflow in {"design", "manifest"} else _important_files(roots)
    for root_index, path, label in matches:
        root = roots[root_index]
        relative = path.relative_to(root).as_posix()
        if workflow in {"design", "manifest"}:
            display_path = _selected_pair_label(root, output_directory)
            file_kind = f"Selected probe pairs · {path.name}"
        else:
            display_path = _display_file_path(path, root, output_directory) if output_directory else relative
            file_kind = f"{label} · {relative}"
        items.append(
            f'<li><a href="{_file_url(token, root_index, root, path)}">{html.escape(display_path)}</a>'
            f'<span class="file-kind">{html.escape(file_kind)}</span></li>'
        )
    return "".join(items) or '<li class="muted">No user-facing output files were found.</li>'


def _primary_map_files(roots: list[Path]) -> list[tuple[int, Path]]:
    maps = _map_files(roots)
    final = [item for item in maps if item[1].name.endswith("_final_probe_map.png")]
    return final or [item for item in maps if item[1].name.endswith("_probe_map.png")]


def _set_job(token: str, **updates: object) -> None:
    with JOBS_LOCK:
        if token in JOBS:
            JOBS[token].update(updates)


def _job_snapshot(token: str) -> dict[str, object] | None:
    with JOBS_LOCK:
        job = JOBS.get(token)
        return dict(job) if job is not None else None


def _active_job_count() -> int:
    with JOBS_LOCK:
        return sum(
            str(job.get("status", "")).lower() in ACTIVE_JOB_STATUSES
            for job in JOBS.values()
        )


def _active_work_exists() -> bool:
    """Return whether queued or running work can still touch shared state."""
    return _active_job_count() > 0 or RUN_LOCK.locked()


def _busy_response(handler: BaseHTTPRequestHandler, action: str) -> None:
    handler.send_response(409)
    body = json.dumps(
        {
            "error": (
                f"Cannot {action} while a design, QC, manifest, plot, or genome-index "
                "run is queued or running. Wait for it to finish or cancel it first."
            )
        }
    ).encode()
    handler.send_header("Content-Type", "application/json")
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def _run_job(
    token: str,
    fields: dict[str, str],
    uploads: list[tuple[str, str, bytes]],
    cancel_event: threading.Event,
) -> None:
    # These are presentation ranges, not scientific estimates. They keep the
    # browser informative during operations whose duration cannot be predicted
    # reliably (network transfer, Bowtie2, and Primer3). Detailed counters,
    # when available, are interpolated within the corresponding phase.
    phase_ranges = {
        "starting": (1, 4),
        "transcript": (4, 12),
        "reference": (12, 18),
        "download": (18, 38),
        "annotation": (38, 48),
        "fasta": (48, 54),
        "premrna": (54, 60),
        "design": (60, 64),
        "designprobes": (64, 78),
        "qc": (78, 88),
        "curation": (88, 94),
        "plot": (94, 97),
        "output": (97, 99),
        "build": (48, 76),
    }
    if fields.get("mode", "").strip().lower() == "index":
        if _field_enabled(fields, "index_annotation_database"):
            phase_ranges["annotation"] = (38, 60)
            phase_ranges["build"] = (60, 99)
        else:
            phase_ranges["build"] = (38, 99)

    last_progress = 1

    def progress(update: dict[str, object]) -> None:
        nonlocal last_progress
        if cancel_event.is_set():
            raise core.RunCancelled("Run cancelled by user")
        phase = str(update.get("phase") or "running")
        completed = update.get("completed")
        total = update.get("total")
        progress_value = None
        phase_start, phase_end = phase_ranges.get(phase, (1, 99))
        if completed is not None and total not in (None, 0):
            try:
                ratio = max(0.0, min(1.0, float(completed) / float(total)))
            except (TypeError, ValueError, ZeroDivisionError):
                ratio = 0.0
            progress_value = round(phase_start + (phase_end - phase_start) * ratio)
        elif phase in phase_ranges:
            progress_value = phase_start
        channel = str(update.get("channel") or "").strip()
        channel_total = update.get("channel_total")
        channel_index = update.get("channel_index")
        if channel and channel_total not in (None, 0) and channel_index not in (None, 0):
            try:
                total_channels = max(1, int(channel_total))
                current_channel = max(1, min(total_channels, int(channel_index)))
                if total_channels > 1:
                    if update.get("channel_fraction") is not None:
                        local_fraction = float(update.get("channel_fraction"))
                    elif progress_value is not None:
                        local_fraction = (float(progress_value) - 1.0) / 98.0
                    else:
                        local_fraction = 0.0
                    local_fraction = max(0.0, min(1.0, local_fraction))
                    progress_value = round(
                        4.0 + 95.0 * ((current_channel - 1 + local_fraction) / total_channels)
                    )
                # A single selected channel is one ordinary design run. Keep
                # its phase progress directly instead of squeezing it through
                # the all-channel aggregation formula.
            except (TypeError, ValueError, ZeroDivisionError):
                progress_value = None
        if progress_value is not None:
            progress_value = max(last_progress, int(progress_value))
            last_progress = progress_value
        message = str(update.get("message") or "Working")
        if channel and channel_total not in (None, 0):
            message = f"{channel} ({channel_index}/{channel_total}) · {message}"
        _set_job(
            token,
            stage=phase,
            message=message,
            progress=progress_value,
        )

    _set_job(token, status="running", stage="starting", message="Starting HCRProbeForge", progress=1)
    try:
        with RUN_LOCK:
            result = _run_submission(fields, uploads, progress_callback=progress, cancel_event=cancel_event)
        success = bool(result.get("success"))
        _set_job(
            token,
            **result,
            status="completed" if success else "failed",
            stage="complete" if success else "failed",
            message="Run completed successfully" if success else "Run finished with an error",
            progress=100,
        )
    except core.RunCancelled:
        _set_job(
            token,
            status="cancelled",
            stage="cancelled",
            message="Run cancelled by user",
            progress=100,
            cancelled=True,
            return_code=None,
            success=False,
            stdout="",
            stderr="Run cancelled by user",
            roots=[],
            output_directory="",
        )
    except Exception as exc:
        _set_job(
            token,
            status="failed",
            stage="failed",
            message=str(exc),
            progress=100,
            return_code=1,
            success=False,
            stdout="",
            stderr=str(exc),
            roots=[],
            output_directory="",
        )


def _progress_page(token: str, request_token: str) -> str:
    token_json = json.dumps(token)
    return f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">{FAVICON_TAG}<title>HCRProbeForge progress</title><style>{STYLE}</style></head><body><main><section class="hero"><div class="brand">{LOGO_SVG}<div><div class="eyebrow">HCR v3 split-initiator probe design</div><h1>Run in progress</h1><div class="subtitle">HCRProbeForge is working. Long designs and multi-gene runs can take several minutes; this page will update automatically.</div></div></div></section><section class="card"><div class="progress-track"><div id="bar" class="progress-bar indeterminate"></div></div><div class="status-line"><strong id="stage">Starting</strong><span id="percent">Working…</span></div><p id="message" class="muted">Preparing the run.</p><p class="small-note">Keep this browser tab open. The output page will appear automatically when the run finishes.</p><form method="post" action="/cancel/{quote(token)}" onsubmit="return confirm('Cancel this run?');"><button class="cancel-button" type="submit">Cancel run</button></form></section></main><script>
const token={token_json};
async function poll(){{try{{const response=await fetch('/status/'+token);const data=await response.json();if(!response.ok)throw new Error(data.error||'Status unavailable');const bar=document.getElementById('bar');const percent=document.getElementById('percent');const progress=data.progress;if(progress===null||progress===undefined){{bar.classList.add('indeterminate');percent.textContent='Working…';}}else{{bar.classList.remove('indeterminate');bar.style.width=progress+'%';percent.textContent=progress+'%';}}document.getElementById('stage').textContent=data.stage||data.status;document.getElementById('message').textContent=data.message||'Working';if(data.status==='completed'||data.status==='failed'||data.status==='cancelled'){{window.hcrNavigating=true;window.location='/result/'+token;return;}}}}catch(error){{document.getElementById('message').textContent=error.message;}}setTimeout(poll,800);}}
poll();</script>{_browser_lifecycle_script(request_token)}</body></html>'''


def _result_page(token: str, job: dict[str, object], request_token: str) -> str:
    roots = _roots_from_job(job)
    success = bool(job.get("success"))
    workflow = str(job.get("workflow", "design"))
    cancelled = job.get("status") == "cancelled" or bool(job.get("cancelled"))
    if cancelled:
        status_text = "cancelled"
        status_detail = "The run was stopped before completion. No new final result was created."
    elif success:
        status_text = "completed"
        status_detail = {
            "design": "Probe design, quality control, and output files were written successfully.",
            "manifest": "The manifest run completed and its results were written successfully.",
            "plot": "The probe map was created successfully.",
            "qc": "The quality-control workbook was written successfully.",
            "index": "The genome reference was downloaded, indexed, and registered successfully.",
        }.get(workflow, "The requested workflow completed successfully.")
    else:
        status_text = "could not complete"
        status_detail = "Review the error details and console output below."
    stderr = str(job.get("stderr", ""))
    error_html = f'<div class="error"><strong>What happened</strong>\n{html.escape(stderr)}</div>' if stderr and not cancelled else ""
    output_directory = str(job.get("output_directory") or (roots[0] if roots else ""))
    maps = _primary_map_files(roots) if workflow in {"design", "manifest", "plot"} else []
    map_cards: list[str] = []
    for root_index, path in maps[:20]:
        root = roots[root_index]
        display_path = _map_display_label(path, root, output_directory)
        map_cards.append(
            f'<figure class="map-card"><img src="{_file_url(token, root_index, root, path)}" alt="{html.escape(path.stem)}" loading="lazy">'
            f'<figcaption>{html.escape(display_path)}</figcaption></figure>'
        )
    map_html = "".join(map_cards) or '<p class="muted">No final probe map was produced for this run.</p>'
    if workflow == "index":
        reference = job.get("reference")
        reference_dict = reference if isinstance(reference, dict) else {}
        index_prefix = str(job.get("index_prefix") or reference_dict.get("index_prefix", ""))
        index_metadata = str(job.get("index_metadata_path") or reference_dict.get("web_metadata_path", ""))
        index_details = "".join(
            part
            for part in (
                f'<p class="small-note">Bowtie2 index prefix</p><code class="path">{html.escape(index_prefix or "Not reported")}</code>',
                f'<p class="small-note">HCRProbeForge index metadata JSON</p><code class="path">{html.escape(index_metadata or "Not reported")}</code>',
            )
        )
        workflow_section = f'<section class="card"><h2>Genome index</h2>{index_details}</section>'
        file_note = "The JSON record stores the organism, assembly, accession, checksum, alias, and registered index location."
    elif workflow in {"design", "manifest", "plot"}:
        workflow_section = f'<section class="card"><h2>Final probe map</h2><div class="map-gallery">{map_html}</div></section>'
        file_note = "Technical and intermediate files are kept under each target folder's <code>details/</code> directory."
    else:
        workflow_section = ""
        file_note = "The QC workbook is written to the folder shown above."
    return f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">{FAVICON_TAG}<title>HCRProbeForge result</title><style>{STYLE}</style></head><body><main><section class="hero"><div class="brand">{LOGO_SVG}<div><div class="eyebrow">HCRProbeForge</div><h1>Run {html.escape(status_text)}</h1><div class="subtitle">{html.escape(status_detail)}</div></div></div></section>{error_html}<section class="card"><div class="result-actions"><a class="button-link" href="/open-directory/{quote(token)}">Open output folder</a><a class="button-link secondary" href="/?new_run={quote(token)}">Start another run</a></div><p class="small-note">Files were written to:</p><code class="path">{html.escape(output_directory or "No output directory")}</code></section>{workflow_section}<section class="card"><details><summary>Console output</summary><div class="console">{html.escape(str(job.get("stdout", "")))}</div></details></section><section class="card"><h2>Key output files</h2><ul class="file-list">{_important_file_links(token, roots, output_directory, workflow)}</ul><p class="small-note">{file_note} Click a file to view or download it; use <strong>Open output folder</strong> to work with the files locally.</p></section></main>{_browser_lifecycle_script(request_token)}</body></html>'''


def _normalise_hostname(value: str | None) -> str:
    """Return a case-insensitive hostname without brackets or a trailing dot."""
    raw = str(value or "").strip()
    if not raw:
        return ""
    try:
        parsed = urlsplit(f"//{raw}")
        hostname = parsed.hostname or ""
    except ValueError:
        hostname = ""
    if not hostname and raw.count(":") > 1:
        hostname = raw.strip("[]")
    return hostname.rstrip(".").lower()


def _is_loopback_bind(host: str) -> bool:
    normalized = _normalise_hostname(host)
    if normalized == "localhost":
        return True
    try:
        return ipaddress.ip_address(normalized).is_loopback
    except ValueError:
        return False


def _allowed_hostnames(bind_host: str, extra_hosts: list[str] | None = None) -> set[str]:
    hosts = {
        _normalise_hostname(bind_host),
        "localhost",
        "127.0.0.1",
        "::1",
        "0.0.0.0",
    }
    hosts.update(_normalise_hostname(value) for value in (extra_hosts or []))
    return {value for value in hosts if value}


class HCRProbeForgeHTTPServer(ThreadingHTTPServer):
    """HTTP server carrying the per-process browser capability settings."""

    daemon_threads = True
    allow_reuse_address = True

    def __init__(
        self,
        server_address: tuple[str, int],
        request_handler: type[BaseHTTPRequestHandler],
        *,
        access_token: str,
        require_auth: bool,
        allowed_hosts: set[str],
    ) -> None:
        super().__init__(server_address, request_handler)
        self.hcr_access_token = access_token
        self.hcr_require_auth = require_auth
        self.hcr_allowed_hosts = allowed_hosts


class HCRProbeForgeHandler(BaseHTTPRequestHandler):
    server_version = f"HCRProbeForge/{core.__version__}"

    @property
    def _hcr_server(self) -> HCRProbeForgeHTTPServer:
        return self.server  # type: ignore[return-value]

    def _send(
        self,
        status: int,
        body: bytes,
        content_type: str = "text/html; charset=utf-8",
        *,
        headers: dict[str, str] | None = None,
    ) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        # Every page and status response is generated from local mutable
        # workflow/reference state. Prevent normal-cache and back/forward
        # cache layers from restoring the pre-index-build form.
        self.send_header("Cache-Control", "no-store, no-cache, must-revalidate, max-age=0")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        for name, value in (headers or {}).items():
            self.send_header(name, value)
        self.end_headers()
        self.wfile.write(body)

    def _cookie_token(self) -> str:
        cookie = SimpleCookie()
        try:
            cookie.load(self.headers.get("Cookie", ""))
        except (ValueError, CookieError):
            return ""
        morsel = cookie.get(WEB_TOKEN_COOKIE)
        return morsel.value.strip() if morsel is not None else ""

    def _request_tokens(self, body_token: str | None = None) -> list[str]:
        values = [self.headers.get(WEB_TOKEN_HEADER, "").strip(), self._cookie_token()]
        if body_token:
            values.append(body_token.strip())
        return [value for value in values if value]

    def _token_matches(self, body_token: str | None = None) -> bool:
        expected = self._hcr_server.hcr_access_token
        return any(secrets.compare_digest(value, expected) for value in self._request_tokens(body_token))

    def _host_is_allowed(self) -> bool:
        raw_host = self.headers.get("Host", "")
        try:
            parsed = urlsplit(f"//{raw_host}")
            port = parsed.port
        except ValueError:
            return False
        hostname = _normalise_hostname(parsed.hostname)
        if not hostname or (hostname not in self._hcr_server.hcr_allowed_hosts and not self._hcr_server.hcr_require_auth):
            return False
        return port == self._hcr_server.server_port

    def _origin_is_allowed(self) -> bool:
        origin = self.headers.get("Origin", "").strip()
        if not origin:
            return True
        if origin.lower() == "null":
            return False
        try:
            parsed = urlsplit(origin)
            port = parsed.port
            request_host = urlsplit(f"//{self.headers.get('Host', '')}")
        except ValueError:
            return False
        hostname = _normalise_hostname(parsed.hostname)
        request_hostname = _normalise_hostname(request_host.hostname)
        expected_port = 443 if parsed.scheme.lower() == "https" else 80 if parsed.scheme.lower() == "http" else None
        return (
            parsed.scheme.lower() in {"http", "https"}
            and hostname == request_hostname
            and (hostname in self._hcr_server.hcr_allowed_hosts or self._hcr_server.hcr_require_auth)
            and port == self._hcr_server.server_port
            and expected_port is not None
        )

    def _auth_cookie_header(self) -> str:
        return f"{WEB_TOKEN_COOKIE}={self._hcr_server.hcr_access_token}; Path=/; SameSite=Strict; HttpOnly"

    def _reject_request(self, status: int = 403, message: str = "Request not authorized") -> bool:
        self._send(status, message.encode(), "text/plain; charset=utf-8")
        return False

    def _authorize_request(self, *, mutation: bool, body_token: str | None = None) -> bool:
        if not self._host_is_allowed():
            return self._reject_request(400, "Host is not allowed")
        if mutation and not self._origin_is_allowed():
            return self._reject_request(403, "Origin is not allowed")
        if (mutation or self._hcr_server.hcr_require_auth) and not self._token_matches(body_token):
            return self._reject_request(403, "Missing or invalid HCRProbeForge access token")
        return True

    def _set_cookie_and_redirect(self, location: str) -> None:
        self.send_response(303)
        self.send_header("Location", location)
        self.send_header("Set-Cookie", self._auth_cookie_header())
        self.send_header("Content-Length", "0")
        self.send_header("Cache-Control", "no-store, no-cache, must-revalidate, max-age=0")
        self.send_header("Pragma", "no-cache")
        self.end_headers()

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if not self._host_is_allowed():
            self._reject_request(400, "Host is not allowed")
            return
        query_token = parse_qs(parsed.query).get("access_token", [""])[0].strip()
        if query_token:
            if parsed.path != "/" or not secrets.compare_digest(query_token, self._hcr_server.hcr_access_token):
                self._reject_request(403, "Invalid HCRProbeForge access token")
                return
            self._set_cookie_and_redirect("/")
            return
        if self._hcr_server.hcr_require_auth and not self._authorize_request(mutation=False):
            return
        if parsed.path == "/favicon.svg":
            self._send(200, FAVICON_SVG.encode(), "image/svg+xml")
            return
        if parsed.path == "/index-status":
            species = parse_qs(parsed.query).get("species", [""])[0].strip()
            preset = references.get_species_preset(species)
            metadata = references.find_installed_reference(species) or {}
            if preset is None and not metadata:
                self._send(404, json.dumps({"error": "Species preset not found"}).encode(), "application/json")
                return
            species_key = preset.key if preset is not None else str(
                metadata.get("index_species_alias") or metadata.get("species") or species
            )
            ready = references.registered_index_is_ready(species_key, metadata=metadata or None)
            self._send(
                200,
                json.dumps(
                    {
                        "species": species_key,
                        "ready": ready,
                        "annotation_ready": _annotation_database_ready(
                            species_key,
                            metadata=metadata or None,
                            # This endpoint feeds the setup-page badge. The
                            # atomic metadata/path state is sufficient for a
                            # fast UI update; the design worker performs the
                            # full SQLite schema validation before querying it.
                            structural_check=False,
                        ),
                        "index_root": str(references.hcrprobedesign_data_root()),
                    }
                ).encode(),
                "application/json",
            )
            return
        if parsed.path == "/delete-species-preview":
            if _active_work_exists():
                _busy_response(self, "inspect or change a species preset")
                return
            species = parse_qs(parsed.query).get("species", [""])[0].strip()
            try:
                plan = references.species_deletion_plan(species)
            except (ValueError, OSError) as exc:
                self._send(400, json.dumps({"error": str(exc)}).encode(), "application/json")
                return
            self._send(200, json.dumps(plan).encode(), "application/json")
            return
        if parsed.path == "/":
            page = _render_form(self._hcr_server.hcr_access_token)
            self._send(200, page.encode(), headers={"Set-Cookie": self._auth_cookie_header()})
            return
        if parsed.path.startswith("/progress/"):
            token = unquote(parsed.path.split("/", 2)[2])
            job = _job_snapshot(token)
            if not job:
                self._send(404, b"Run not found", "text/plain")
                return
            if job.get("status") in {"completed", "failed", "cancelled"}:
                self.send_response(303)
                self.send_header("Location", f"/result/{quote(token)}")
                self.end_headers()
                return
            self._send(200, _progress_page(token, self._hcr_server.hcr_access_token).encode())
            return
        if parsed.path.startswith("/status/"):
            token = unquote(parsed.path.split("/", 2)[2])
            job = _job_snapshot(token)
            if not job:
                self._send(404, json.dumps({"error": "Run not found"}).encode(), "application/json")
                return
            status = {
                "status": job.get("status", "queued"),
                "stage": job.get("stage", "queued"),
                "message": job.get("message", "Queued"),
                "progress": job.get("progress"),
                "cancel_requested": job.get("cancel_requested", False),
            }
            self._send(200, json.dumps(status).encode(), "application/json")
            return
        if parsed.path.startswith("/result/"):
            token = unquote(parsed.path.split("/", 2)[2])
            job = _job_snapshot(token)
            if not job:
                self._send(404, b"Result not found", "text/plain")
                return
            if job.get("status") not in {"completed", "failed", "cancelled"}:
                self.send_response(303)
                self.send_header("Location", f"/progress/{quote(token)}")
                self.end_headers()
                return
            self._send(200, _result_page(token, job, self._hcr_server.hcr_access_token).encode())
            return
        if parsed.path.startswith("/open-directory/"):
            token = unquote(parsed.path.split("/", 2)[2])
            job = _job_snapshot(token)
            if not job:
                self._send(404, b"Run not found", "text/plain")
                return
            roots = _roots_from_job(job)
            directory_value = str(job.get("output_directory") or "").strip()
            directory = Path(directory_value).expanduser() if directory_value else Path()
            if not directory.is_dir() and roots:
                directory = roots[0]
            if not directory_value and not roots:
                self._send(404, b"No output directory is available", "text/plain")
                return
            try:
                _open_local_directory(directory)
            except (FileNotFoundError, RuntimeError) as exc:
                self._send(409, str(exc).encode(), "text/plain")
                return
            self.send_response(303)
            self.send_header("Location", f"/result/{quote(token)}")
            self.end_headers()
            return
        if parsed.path.startswith("/file/"):
            pieces = parsed.path.split("/", 3)
            if len(pieces) != 4:
                self._send(404, b"File not found", "text/plain")
                return
            token = unquote(pieces[2])
            job = _job_snapshot(token)
            if not job:
                self._send(404, b"File not found", "text/plain")
                return
            tail = pieces[3].split("/", 1)
            if len(tail) != 2:
                self._send(404, b"File not found", "text/plain")
                return
            try:
                root = _roots_from_job(job)[int(tail[0])]
            except (IndexError, ValueError):
                self._send(404, b"File not found", "text/plain")
                return
            base = root.resolve()
            target = (base / unquote(tail[1])).resolve()
            if not target.is_relative_to(base) or not target.is_file():
                self._send(403, b"Invalid file", "text/plain")
                return
            self._send(200, target.read_bytes(), mimetypes.guess_type(target.name)[0] or "application/octet-stream")
            return
        self._send(404, b"Not found", "text/plain")

    def do_POST(self) -> None:
        endpoint = urlparse(self.path).path
        body_token = None
        if endpoint == "/shutdown":
            try:
                body_length = int(self.headers.get("Content-Length", "0"))
            except ValueError:
                self._reject_request(400, "Invalid Content-Length")
                return
            if body_length > 4096:
                self._reject_request(413, "Shutdown request is too large")
                return
            if body_length:
                body_token = self.rfile.read(body_length).decode("utf-8", errors="replace").strip()
        if not self._authorize_request(mutation=True, body_token=body_token):
            return
        if endpoint == "/shutdown":
            with LIFECYCLE_LOCK:
                if _active_work_exists():
                    _busy_response(self, "shut down the local webapp")
                    return
                self._send(204, b"", "text/plain; charset=utf-8")
            threading.Thread(target=self.server.shutdown, daemon=True).start()
            return
        if endpoint == "/open-index-folder":
            try:
                index_root = references.hcrprobedesign_data_root()
                index_root.mkdir(parents=True, exist_ok=True)
                _open_local_directory(index_root)
            except (FileNotFoundError, OSError, RuntimeError) as exc:
                self._send(409, json.dumps({"error": str(exc)}).encode(), "application/json")
                return
            self._send(200, json.dumps({"path": str(index_root)}).encode(), "application/json")
            return
        if endpoint == "/choose-folder":
            selected = _choose_local_folder()
            self._send(200, json.dumps({"path": selected}).encode(), "application/json")
            return
        if endpoint.startswith("/cancel/"):
            token = unquote(endpoint.split("/", 2)[2])
            job = _job_snapshot(token)
            if not job:
                self._send(404, b"Run not found", "text/plain")
                return
            if job.get("status") in {"completed", "failed", "cancelled"}:
                self.send_response(303)
                self.send_header("Location", "/")
                self.end_headers()
                return
            cancel_event = job.get("cancel_event")
            if isinstance(cancel_event, threading.Event):
                cancel_event.set()
            _set_job(token, stage="cancelling", message="Cancelling run…", cancel_requested=True)
            self.send_response(303)
            self.send_header("Location", "/")
            self.end_headers()
            return
        try:
            fields, uploads = _parse_multipart(self)
            if endpoint == "/register-species":
                with LIFECYCLE_LOCK:
                    if _active_work_exists():
                        _busy_response(self, "add a species preset")
                        return
                    result = _register_species_submission(fields, uploads)
                self._send(200, json.dumps(result).encode(), "application/json")
                return
            if endpoint == "/delete-species":
                with LIFECYCLE_LOCK:
                    if _active_work_exists():
                        _busy_response(self, "remove a species preset")
                        return
                    result = _delete_species_submission(fields)
                self._send(200, json.dumps(result).encode(), "application/json")
                return
            if endpoint == "/inspect":
                self._send(200, json.dumps(_inspect_submission(fields, uploads)).encode(), "application/json")
                return
            if endpoint != "/run":
                self._send(404, b"Not found", "text/plain")
                return
            token = uuid4().hex
            cancel_event = threading.Event()
            with LIFECYCLE_LOCK:
                with JOBS_LOCK:
                    JOBS[token] = {
                        "workflow": fields.get("mode", "design").strip().lower() or "design",
                        "status": "queued",
                        "stage": "queued",
                        "message": "Waiting for the local worker",
                        "progress": 0,
                        "return_code": None,
                        "success": False,
                        "stdout": "",
                        "stderr": "",
                        "roots": [],
                        "output_directory": "",
                        "cancel_event": cancel_event,
                        "cancel_requested": False,
                        "cancelled": False,
                    }
            worker = threading.Thread(target=_run_job, args=(token, fields, uploads, cancel_event), daemon=True)
            worker.start()
            self.send_response(303)
            self.send_header("Location", f"/progress/{quote(token)}")
            self.end_headers()
        except Exception as exc:
            if endpoint in {"/register-species", "/delete-species"}:
                self._send(400, json.dumps({"error": str(exc)}).encode(), "application/json")
                return
            page = f'<!doctype html><html><head><meta charset="utf-8">{FAVICON_TAG}<style>{STYLE}</style></head><body><main><h1>Run could not start</h1><div class="error">{html.escape(str(exc))}</div><p><a href="/">Return to setup</a></p></main></body></html>'
            self._send(400, page.encode())

    def log_message(self, format: str, *args: object) -> None:
        if os.environ.get("HCRPROBEFORGE_WEB_VERBOSE"):
            super().log_message(format, *args)


def _port_value(value: str) -> int:
    try:
        port = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("port must be a whole number") from exc
    if not 0 <= port <= 65535:
        raise argparse.ArgumentTypeError("port must be between 0 and 65535")
    return port


def _validate_auth_token(value: str) -> str:
    token = value.strip()
    if len(token) < 16 or any(character.isspace() or character in ";,\r\n" for character in token):
        raise ValueError("--auth-token must contain at least 16 non-whitespace characters")
    return token


def _bind_web_server(
    host: str,
    requested_port: int,
    *,
    access_token: str,
    require_auth: bool,
    allowed_hosts: set[str],
) -> HCRProbeForgeHTTPServer:
    """Bind the requested port, trying nearby ports when it is occupied."""
    if requested_port == 0:
        ports = [0]
    else:
        ports = list(range(requested_port, min(65535, requested_port + WEB_PORT_ATTEMPTS - 1) + 1))
    last_error: OSError | None = None
    for port in ports:
        try:
            return HCRProbeForgeHTTPServer(
                (host, port),
                HCRProbeForgeHandler,
                access_token=access_token,
                require_auth=require_auth,
                allowed_hosts=allowed_hosts,
            )
        except OSError as exc:
            if exc.errno != errno.EADDRINUSE:
                raise
            last_error = exc
    tried = ", ".join(str(port) for port in ports)
    raise OSError(errno.EADDRINUSE, f"Could not bind HCRProbeForge on ports {tried}") from last_error


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Launch the local HCRProbeForge webapp")
    parser.add_argument("--host", default="127.0.0.1", help="Interface to bind; default is local-only")
    parser.add_argument("--port", type=_port_value, default=DEFAULT_WEB_PORT)
    parser.add_argument(
        "--auth-token",
        help="Access token for non-loopback binding; remote mode requires at least 16 characters",
    )
    parser.add_argument(
        "--allowed-host",
        action="append",
        default=[],
        metavar="HOST",
        help="Additional Host header/name allowed for the webapp; may be repeated",
    )
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args(argv)
    require_auth = not _is_loopback_bind(args.host)
    if require_auth and not args.auth_token:
        parser.error("--auth-token is required when --host is not loopback")
    try:
        access_token = _validate_auth_token(args.auth_token) if args.auth_token else secrets.token_urlsafe(32)
    except ValueError as exc:
        parser.error(str(exc))
    allowed_hosts = _allowed_hostnames(args.host, args.allowed_host)
    try:
        server = _bind_web_server(
            args.host,
            args.port,
            access_token=access_token,
            require_auth=require_auth,
            allowed_hosts=allowed_hosts,
        )
    except OSError as exc:
        parser.error(str(exc))
    url = f"http://{args.host}:{server.server_port}/"
    if server.server_port != args.port and args.port != 0:
        print(f"Port {args.port} is unavailable; using port {server.server_port}.")
    if require_auth:
        access_url = f"{url}?access_token={quote(access_token, safe='')}"
        print(f"HCRProbeForge webapp is running at {url}")
        print(f"Access URL: {access_url}")
    else:
        print(f"HCRProbeForge webapp is running at {url}")
    if not args.no_browser and (not require_auth or args.host not in {"0.0.0.0", "::"}):
        browser_url = access_url if require_auth else url
        threading.Timer(0.6, lambda: _open_browser_url(browser_url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
