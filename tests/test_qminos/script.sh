#!/usr/bin/bash
set -e

activate-python3-venv ()
{
    source /opt/python3-venv/$1/bin/activate
}

activate-python3-venv p37-coralME && python3 script.py 1> p37.out 2> p37.err && deactivate
activate-python3-venv p38-coralME && python3 script.py 1> p38.out 2> p38.err && deactivate
activate-python3-venv p39-coralME && python3 script.py 1> p39.out 2> p39.err && deactivate
activate-python3-venv p310-coralME && python3 script.py 1> p310.out 2> p310.err && deactivate
activate-python3-venv p311-coralME && python3 script.py 1> p311.out 2> p311.err && deactivate
activate-python3-venv p312-coralME && python3 script.py 1> p312.out 2> p312.err && deactivate
activate-python3-venv p313-coralME && python3 script.py 1> p313.out 2> p313.err && deactivate
