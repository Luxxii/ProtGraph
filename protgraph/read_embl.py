import csv

import io
import mmap

from Bio import SwissProt

def rows(f, chunksize=4096):
    """ 
    Split entry and add as bytearray into the queue 
    """
    curr_row = bytearray()
    sep_array = bytearray(b"\nID")
    beg_array = bytearray(b"ID")

    # Alternative, which only uses read/find and tell
    # This is slower but reads the file faster
    # sep=b'\nID'
    # while True:
    #     pos = f.find(sep)
    #     if pos == -1:
    #         yield f.read()
    #         break
    #     yield f.read(pos - f.tell() + 1)

    # The Graph Generator parse the entries
    # offset = 0
    for chunk in iter(lambda: f.read(chunksize), b''):
        curr_row.extend(chunk)
        while True:
            i = curr_row.split(sep_array, 1)  # If python would allow us to add an offset, this could be even faster!
            if len(i) == 1:
                break  # We could update here the offset
            yield i[0]
            curr_row = beg_array + i[1]


def read_embl(path_to_embls: list, queue):
    """ Reads entries from a list of existing embl files """
    for input_f in path_to_embls:
        with open(input_f, "rb") as in_file:
            with mmap.mmap(in_file.fileno(), 0, access=mmap.ACCESS_READ) as mf:
                for r in rows(mf): 
                    queue.put(r)
