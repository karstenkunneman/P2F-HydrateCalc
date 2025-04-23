import streamlit

import streamlit.web.cli as stcli
import os, sys

import numpy
import streamlit
import streamlit.runtime.scriptrunner.magic_funcs
import matplotlib
import thermo
import openpyxl

def resolve_path(path):
    resolved_path = os.path.abspath(os.path.join(os.getcwd(), path))
    return resolved_path


if __name__ == "__main__":
    sys.argv = [
        "streamlit",
        "run",
        resolve_path("app.py"),
        "--global.developmentMode=false",
    ]
    sys.exit(stcli.main())