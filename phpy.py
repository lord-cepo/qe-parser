from utils.pwpy_utils import XmlQe, py2fort

def py2ph(py_dict, outfile='py.ph.in', first_line='generated with py2ph'):
    lines = [f'{first_line}\n', '&INPUTPH\n']
    for key, value in py_dict.items():
        t = XmlQe('PH')
        var = t.search_by_name(key)
        tipo = var.get('type') 
        if tipo is None:
            tipo = t.parent_map[var].get('type')
        if not isinstance(value, XmlQe.conversions[tipo]):
            raise ValueError(f'Variable {key} should be of type {var.get("type")}')
        options = [v2.strip().replace("'", '') for v in var.findall('.//opt') for v2 in v.get('val').split(',')]
        if len(options) > 0 and value not in options:
            raise ValueError(f'Variable {value} should be one of {options}')
        lines.append(f'  {key} = {py2fort(value)}\n')
    lines.append('/\n\n')
    with open(outfile, 'w', encoding='utf-8') as f:
        f.writelines(lines)

    