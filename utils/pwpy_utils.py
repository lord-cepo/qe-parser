from typing import Union
import os

def py2fort(data_value):
    if type(data_value) is bool:
        return '.true.' if data_value else '.false.'
    if type(data_value) is int:
        return str(data_value)
    if type(data_value) is float:
        return str(data_value).replace('e', 'd')
    if type(data_value) is str:
        return f"'{data_value}'"
    if type(data_value) is list:
        return ' '.join(map(str, data_value))

def fort2py(data_type, data_value):
    if data_type == 'INTEGER':
        return int(data_value)
    if data_type == 'REAL':
        return float(data_value.replace('d', 'e'))
    if data_type == 'CHARACTER':
        return data_value.replace("'", "")
    if data_type == 'LOGICAL':
        return data_value == '.true.'
    raise ValueError(f'Unknown data type {data_type}')


class XmlQe:
    conversions = {
        'INTEGER': int,
        'REAL': Union[float, int],
        'CHARACTER': str,
        'LOGICAL': bool,
    }
    mass = {
        'Bi': 208.9804,
        'Te': 127.6,
    }
    def __init__(self, exe, filename=None):
        import xml.etree.ElementTree as ET
        if filename is None:
            filename = os.path.dirname(os.path.abspath(__file__))+'/INPUT_'+exe+'.xml'
        self.tree = ET.parse(filename)
        self.root = self.tree.getroot()
        self.namelists = {n.get('name'):[] for n in self.root.findall('.//namelist')}
        self.parent_map = {c: p for p in self.tree.iter() for c in p}
    
    def search_by_name(self, name):
        if "(" in name:
            xml_type = 'dimension'
        else:
            xml_type = 'var'
        key_to_search = name.split("(")[0]
        node = self.root.find(f'.//{xml_type}[@name="{key_to_search}"]')
        if node is None:
            raise ValueError(f'Variable {key_to_search} not found')
        return node