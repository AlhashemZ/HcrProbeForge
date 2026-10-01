from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from hcrprobeforge import webapp


class PlatformIntegrationTests(unittest.TestCase):
    def test_wsl_browser_uses_windows_explorer(self):
        with patch.object(webapp, "_is_wsl", return_value=True), patch.object(
            webapp.shutil, "which", return_value="/mnt/c/Windows/explorer.exe"
        ), patch.object(webapp.subprocess, "Popen") as popen, patch.object(
            webapp.webbrowser, "open", return_value=False
        ):
            self.assertTrue(webapp._open_browser_url("http://127.0.0.1:8766/"))
            popen.assert_called_once_with(
                ["/mnt/c/Windows/explorer.exe", "http://127.0.0.1:8766/"],
                stdout=webapp.subprocess.DEVNULL,
                stderr=webapp.subprocess.DEVNULL,
                start_new_session=True,
            )

    def test_wsl_folder_button_converts_path_for_windows_explorer(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            with patch.object(webapp, "_is_wsl", return_value=True), patch.object(
                webapp.shutil, "which", side_effect=lambda name: {
                    "explorer.exe": "/mnt/c/Windows/explorer.exe",
                    "wslpath": "/usr/bin/wslpath",
                }.get(name)
            ), patch.object(
                webapp.subprocess, "run", return_value=subprocess.CompletedProcess(
                    ["wslpath"], 0, "C:\\Users\\tester\\hcr_results\n", ""
                ),
            ), patch.object(webapp.subprocess, "Popen") as popen:
                webapp._open_local_directory(directory)
                popen.assert_called_once_with(
                    ["/mnt/c/Windows/explorer.exe", "C:\\Users\\tester\\hcr_results"],
                    stdout=webapp.subprocess.DEVNULL,
                    stderr=webapp.subprocess.DEVNULL,
                    start_new_session=True,
                )

    def test_non_wsl_browser_keeps_standard_browser_behavior(self):
        with patch.object(webapp, "_is_wsl", return_value=False), patch.object(
            webapp.webbrowser, "open", return_value=True
        ) as browser:
            self.assertTrue(webapp._open_browser_url("http://127.0.0.1:8766/"))
            browser.assert_called_once_with("http://127.0.0.1:8766/")


if __name__ == "__main__":
    unittest.main()
