import subprocess
import re


def get_snps(line):
    if line=='NONE':
        return []
    else:
        full_list = list(map(lambda x: re.sub('\([\d]+\)', '', x), line.split(',')))
        chr_bp_list = list(map(get_ch_bp, full_list))
        return chr_bp_list
    
    

def run(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
    output, error = process.communicate()
    return output, error



def flatten(l):
    return [item for sublist in l for item in sublist]




def get_ch_bp(x):
    return ':'.join(x.split(":")[:2])