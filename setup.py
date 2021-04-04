"""
setup.py

Arrange simulator into a single .exe file
"""

'''----------------------------------------------------------------------
                                IMPORTS
----------------------------------------------------------------------'''

import sys
from cx_Freeze import setup, Executable


'''----------------------------------------------------------------------
                                SETTINGS
----------------------------------------------------------------------'''

base = "Console"
if sys.platform == "win32":
    base = "Win32GUI"

options = {
    "build_exe": {
        # Sometimes a little fine-tuning is needed
        # exclude all backends except wx
        "excludes": ["gtk", "PyQt4", "PyQt5", "tkinter","pandas","tcl","sympy"]
    }
}

executables = [
    Executable("rims.py")
]

setup(
    name="rims0221",
    version="0.1",
    description="ratchet based ion movement simulator",
    executables=executables,
    options=options,
)