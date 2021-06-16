""" hdf5_Partition.py
Usage:
    hdf5_Partition.py <Index> <fileIn> <fileOut>

"""
import numpy as np
import h5py as h5
import pathlib
import sys
#########################################################################################
"""
    The purpose of this script is to copy and partition hdf5 files in the dedalus format
"""
#########################################################################################
def descend_obj(obj,sep = '\t'):
    """
    Iterate through groups in a HDF5 file and prints the groups and datasets names and datasets attributes
    """
    if type(obj) in [h5._hl.group.Group, h5._hl.files.File]:
        for key in obj.keys():
            print(sep, '-', key, ':', obj[key])
            descend_obj(obj[key], sep = sep+'\t')
    elif type(obj) == h5._hl.dataset.Dataset:
        pass
        """for key in obj.attrs.keys():
            print(sep + '\t', '-', key, ':', obj.attrs[key])"""

def copy_attributes(in_object, out_object):
    '''Copy attributes between 2 HDF5 objects.'''
    for key, value in in_object.attrs.items():
        out_object.attrs[key] = value
        print(value)
    return None

def dedalus_parent_attributes(index, in_object, out_object, handler_name = None):
    '''Copy attributes between 2 HDF5 objects.'''
    for key, value in in_object.attrs.items():
        if (key == 'handler_name'):
            print(value)
            if (handler_name == None):
                out_object.attrs[key] = value + '_parition'
            else:
                out_object.attrs[key] = str(handler_name)
            print(out_object.attrs[key])
        elif (key == 'writes'):
            print(value)
            print(index)
            out_object.attrs[key] = index
        else:
            out_object.attrs[key] = value
    return None

def dataset_Group_Name(full_Path, key_Value):
    """
        Output the group path of key_Value
    """
    length = len(key_Value)
    group_Path = full_Path[:-length]
    
    return group_Path

def update_Tuple(tuple_Object, value_Update):
    ################################################################################
    """
        Function that takes a tuple and changes values based on value_Update;
        :::
            'None' corresponds to no change
    """
    ################################################################################
    
    change_Tuple = list(tuple_Object)
    try:
       len(tuple_Object) == len(value_Update)
    except:
        print('Check the lengths of tuple_Object and value_Update tuple')
    for i, value in enumerate(value_Update):
        if (value == None):
            pass
        else:
            change_Tuple[i] = value
    change_Tuple = tuple(change_Tuple)
    
    return change_Tuple
    
def partition_Structure(index, h5file1, h5file2):
    '''Goes through entire h5file1 structure and copy into h5file2 with index alteration.'''
    dedalus_parent_attributes(index, h5file1, h5file2, handler_name = handler_name)
    for key in h5file1.keys():
        if (type(h5file1[key]) == h5._hl.group.Group):
            group_path = h5file1[key].name    
            h5file2.require_group(group_path)
            copy_attributes(h5file1[key], h5file2[key])
            partition_Structure(index, h5file1[key], h5file2[key])

        elif (type(h5file1[key]) == h5._hl.dataset.Dataset):
            name = h5file1[key].name
            group_Path = dataset_Group_Name(name, key)
            dtype = h5file1[key].dtype
            shape = h5file1[key].shape
            maxshape = h5file1[key].maxshape
            chunks = h5file1[key].chunks
            if (len(shape) == 1):
                if (maxshape[0] != None):
                    data = h5file1[key][()]
                    h5file2[group_Path].create_dataset(name = key, data = data, dtype = dtype, 
                                                shape = shape, maxshape = maxshape, chunks = chunks)
                    copy_attributes(h5file1[key], h5file2[key])
                else:
                    shape = update_Tuple(shape, value_Update = (index,))
                    data = h5file1[key][()][:index]
                    h5file2[group_Path].create_dataset(name = key, data = data, dtype = dtype, 
                                                shape = shape, maxshape = maxshape, chunks = chunks)
                    copy_attributes(h5file1[key], h5file2[key])
                    
            elif (len(shape) == 4):
                shape = update_Tuple(shape, value_Update = (index, None, None, None))
                maxshape = update_Tuple(shape, value_Update = (index, None, None, None))
                data = h5file1[key][()][:index, :, :, :]
                h5file2[group_Path].create_dataset(name = key, data = data, dtype = dtype, 
                                            shape = shape, maxshape = maxshape, chunks = chunks)
                copy_attributes(h5file1[key], h5file2[key])
    return None
            
        
        
            
parent_file = pathlib.Path(sys.argv[-2])     
child_file = pathlib.Path(sys.argv[-1])
index = int(sys.argv[-3])
pf = h5.File(parent_file, 'r')
cf = h5.File(child_file, 'w')

partition_Structure(index, pf, cf)

pf.close()
cf.close()
