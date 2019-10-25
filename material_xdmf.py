# SUMMARY:      material_xdmf.py
# USAGE:        collect unweighted rtd
# ORG:          Pacific Northwest National Laboratory
# AUTHOR:       Xuehang Song
# E-MAIL:       xuehang.song@pnnl.gov
# ORIG-DATE:    July-2018
# DESCRIPTION:  revised from "pt_residence_time_flux_v4_hanford.py"
# DESCRIPTION-END

from xml.etree import ElementTree as ET
from xml.dom import minidom
import argparse
import numpy as np
import h5py as h5

# format XDMF outputs


def prettify(element):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ET.tostring(element, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")


def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-dir',
                        '--output_dir',
                        type=str,
                        default="/global/cscratch1/sd/xhsong/npt/base/pt/0/")
    parser.add_argument('-loc_range',
                        '--loc_range',
                        type=str,
                        default="0-9999")
    parser.add_argument('-time_range',
                        '--time_range',
                        type=str,
                        default="160-180")  # 140~160
    parser.add_argument('-dt',
                        '--dt',
                        type=float,
                        default=1.)

    args = vars(parser.parse_args())
    return(args)


coordinate_hdf5 = "/media/sf_e/john/optim_5/single_h5/Coordinates.h5"
material_hdf5 = "/media/sf_e/john/optim_5/300A_material.h5"
output_dir = "/media/sf_e/john/particles/mapped/"

data_hdf5 = "material_data.h5"
materials_xdmf = "materials.xdmf"


# read coordinates
hdf5 = h5.File(coordinate_hdf5, "r")
x = np.array(hdf5["Coordinates"]['X [m]'])+593000
y = np.array(hdf5["Coordinates"]['Y [m]'])+114500
z = np.array(hdf5["Coordinates"]['Z [m]'])
ox = x[0]
oy = y[0]
oz = z[0]
ex = x[-1]
ey = y[-1]
ez = z[-1]
dx = np.diff(x)
dy = np.diff(y)
dz = np.diff(z)
nx = len(dx)
ny = len(dy)
nz = len(dz)
x = ox+np.cumsum(dx)-0.5*dx
y = oy+np.cumsum(dy)-0.5*dy
z = oz+np.cumsum(dz)-0.5*dz
hdf5.close()

# read materials
hdf5 = h5.File(material_hdf5, "r")
material_ids = np.array(hdf5['Materials']["Material Ids"])
material_ids = material_ids.reshape((nx, ny, nz), order="C")
hdf5.close()


hdf5 = h5.File(output_dir+data_hdf5, "w")
hdf5.create_dataset("Materials", data=material_ids)
hdf5.close()

xml_root = ET.Element("Xdmf", Version="3.0")
xml_domain = ET.SubElement(xml_root, "Domain")


xml_grid = ET.SubElement(xml_domain, "Grid",
                         {'Name': "Material",
                          'GridType': 'Uniform'})

xml_toplogoy = ET.SubElement(xml_grid, "Topology",
                             {'TopologyType': '3DRECTMesh',
                              'Dimensions': "{0} {1} {2}".format(nz, ny, nx)})
xml_geometry = ET.SubElement(xml_grid, 'Geometry',
                             {'GeometryType': "VXVYVZ"})
xml_geometry_x = ET.SubElement(xml_geometry, 'DataItem',
                               {'Dimensions': str(nx),
                                "NumberType": "Float",
                                "Precision": "8",
                                "Format": "XML"})
xml_geometry_x.text = np.array_str(x).strip("[]").replace("\n", " ")
xml_geometry_y = ET.SubElement(xml_geometry, 'DataItem',
                               {'Dimensions': str(ny),
                                "NumberType": "Float",
                                "Precision": "8",
                                "Format": "XML"})
xml_geometry_y.text = np.array_str(y).strip("[]").replace("\n", " ")
xml_geometry_z = ET.SubElement(xml_geometry, 'DataItem',
                               {'Dimensions': str(nz),
                                "NumberType": "Float",
                                "Precision": "8",
                                "Format": "XML"})
xml_geometry_z.text = np.array_str(z).strip("[]").replace("\n", " ")

xml_material = ET.SubElement(xml_grid, "Attribute",
                             {"Name": "Materials",
                              "AttributeType": "Scalar",
                              "Center": "Node"})
material_dataitem = ET.SubElement(xml_material, "DataItem",
                                  {"Format": "HDF",
                                   "NumberType": "Float",
                                   "Precision": "8",
                                   "Dimensions": "{0} {1} {2}".format(nx, ny, nz)})
material_dataitem.text = data_hdf5+":/Materials"
fname = output_dir+materials_xdmf
with open(fname, 'w') as f:
    f.write(prettify(xml_root))
