import sys
from cx_Freeze import setup, Executable

"""
setup.py

Arrange simulator into a single .exe file
"""

'''----------------------------------------------------------------------
                                SETTINGS
----------------------------------------------------------------------'''

base = "Console"
if sys.platform == "win32":
    base = "Win32GUI"

options = {
    "build_exe": {
        # exclude all backends except numpy & matplotlib
        "excludes": ["gtk", "PyQt4", "PyQt5", "tkinter", "pandas", "tcl", "sympy", "pgk_resources"]
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
