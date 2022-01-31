import numpy as np
import optparse
from os import listdir
from os.path import isfile, join
import json
import yaml

def read_inputs(input_dir, fw_config):
    stimuli = {}
    input_files = [f for f in listdir(input_dir) if isfile(join(input_dir,f))]
    nbx = -1
    for f in input_files:
        modhash = "".join((filter(str.isdigit, f)))
        save_module = True
        if fw_config:
            with open(fw_config, 'r') as cfg_file:
                params = yaml.safe_load(cfg_file)
            if 'modules' in params:
                save_module = (modhash in params['modules'])
        if save_module:
            stimuli[modhash] = np.loadtxt(
                join(input_dir,f),
                dtype = 'str',
                usecols = range(0,48) # always 48 TCs per module
            )
            if nbx == -1:
                nbx = len(stimuli[modhash])
    return stimuli, nbx


def order_per_bx(stimuli, Nbx):
    stimuli_per_bx = {'event' : []}
    for bx in range(Nbx):
        stimuli_per_bx['event'].append({})
        stimuli_per_bx['event'][bx]['module'] = {}
        for modhash in stimuli:
            stimuli_per_bx['event'][bx]['module'][modhash] = {}
            stimuli_per_bx['event'][bx]['module'][modhash]['tc'] = []
            for tc_e in stimuli[modhash][bx,:]:
                stimuli_per_bx['event'][bx]['module'][modhash]['tc'].append(tc_e)
    return stimuli_per_bx


def main(opt):
    mod_stimuli, nbx = read_inputs(opt.input_dir, opt.fw_config)
    mod_stimuli = order_per_bx(mod_stimuli, nbx)
    output_file = opt.output_dir+'/TClist.json'
    with open(output_file,'w') as f_out:
        json.dump(mod_stimuli, f_out, indent=4)

if __name__=='__main__':

    parser = optparse.OptionParser()
    parser.add_option("--fwcfg", type="string", dest="fw_config",  default="", help="(optional) yaml file with TC->bin mapping information")
    parser.add_option("--in",  type="string", dest="input_dir", help="Input directory")
    parser.add_option("--out", type="string", dest="output_dir", help="Output directory")
    opt, args = parser.parse_args()

    main(opt)
