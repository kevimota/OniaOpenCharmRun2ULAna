import re
import awkward as ak

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def build_p4(acc):
    p4 = ak.zip({'x': acc['x'].value, 
                 'y': acc['y'].value,
                 'z': acc['z'].value,
                 't': acc['t'].value}, with_name="LorentzVector")

    return p4

def build_acc(obj, name=None):
    from coffea import processor
    acc = processor.dict_accumulator({})
    for var in obj.fields:
        acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(obj[var])))
    if name is not None:
        acc[f"n{name}"] = processor.column_accumulator(ak.to_numpy(ak.num(obj)))
    return acc