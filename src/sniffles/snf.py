#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created: 26.09.2021
# Author:  Moritz Smolka
# Contact: moritz.g.smolka@gmail.com
#

import pickle
import json
import gzip
import math

from sniffles import sv

class SNFile:
    def __init__(self,config,handle,filename=None):
        self.config=config
        self.handle=handle
        self.filename=filename
        self.blocks={}
        self.index={}
        self.total_length=0

    def store(self,svcand):
        block_index=int(svcand.pos/self.config.snf_block_size)*self.config.snf_block_size
        if not block_index in self.blocks:
            self.blocks[block_index]={svtype: [] for svtype in sv.TYPES}
            self.blocks[block_index]["_COVERAGE"]={}
        if not self.config.output_rnames:
            svcand.rnames=None
        self.blocks[block_index][svcand.svtype].append(svcand)

    def annotate_block_coverages(self,lead_provider,resolution=500):
        config=self.config
        start_bin=lead_provider.covrtab_min_bin
        end_bin=int(lead_provider.end/config.coverage_binsize)*config.coverage_binsize
        coverage_fwd=0
        coverage_rev=0

        coverage_sum=0
        bin_count=0

        coverage_binsize_combine=self.config.coverage_binsize_combine
        snf_block_size=config.snf_block_size

        for bin in range(start_bin, end_bin+config.coverage_binsize,config.coverage_binsize):
            if bin in lead_provider.covrtab_fwd:
                coverage_fwd+=lead_provider.covrtab_fwd[bin]

            if bin in lead_provider.covrtab_rev:
                coverage_rev+=lead_provider.covrtab_rev[bin]

            coverage_sum+=coverage_fwd+coverage_rev
            bin_count+=1

            if bin%coverage_binsize_combine==0:
                block_index=int(bin/snf_block_size)*snf_block_size

                coverage_total_curr=math.ceil(coverage_sum/float(bin_count))
                if coverage_total_curr > 0:
                    if not block_index in self.blocks:
                        self.blocks[block_index]={svtype: [] for svtype in sv.TYPES}
                        self.blocks[block_index]["_COVERAGE"]={}

                    self.blocks[block_index]["_COVERAGE"][bin]=coverage_total_curr

                coverage_sum=0
                bin_count=0

    def serialize_block(self,block_id):
        return pickle.dumps(self.blocks[block_id])

    def unserialize_block(self,data):
        return pickle.loads(data)

    def write_and_index(self):
        offset=0
        for block_id in sorted(self.blocks):
            data=gzip.compress(self.serialize_block(block_id))
            self.handle.write(data)
            data_len=len(data)
            self.index[block_id]=(offset,data_len)
            offset+=data_len
            self.total_length+=data_len

    def read_header(self):
        try:
            header_text=self.handle.readline()
            self.header_length=len(header_text)
            self.header=json.loads(header_text.strip())
        except Exception as e:
            print(f"Error when reading SNF header from '{self.filename}': {e}. The file may not be a valid .snf file or could have been corrupted.")
            raise e
        self.index=self.header["index"]

    def read_blocks(self,contig,block_index):
        block_index=str(block_index)
        if not contig in self.index:
            return None

        if not block_index in self.index[contig]:
            return None

        blocks=[]
        for block_data_start,block_data_length in self.index[contig][block_index]:
            try:
                self.handle.seek(self.header_length+block_data_start)
                data=gzip.decompress(self.handle.read(block_data_length))
                blocks.append(self.unserialize_block(data))
            except Exception as e:
                print(f"Error when reading block '{contig}.{block_index}' from '{self.filename}': {e}. The file may not be a valid .snf file or could have been corrupted.")
                raise e
        return blocks

    def get_index(self):
        return self.index

    def get_total_length(self):
        return self.total_length

    def close(self):
        self.handle.close()
