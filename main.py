from pwpy import pw2py, py2pw

pydict = pw2py('examples/pw.in')

py2pw(pydict, 'examples/py.pw.in')

