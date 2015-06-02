#! read and get information of gff file
#  made by SSW

import re

class gff_items_class():
    def __init__(self,line,types_included=[],target_regexps=['ID']):
        self.line=line.rstrip('\r\n')
        self.types_included=types_included
        self.target_regexps=target_regexps
    def parse_gff_line(self):
        items_dict = {}
        line_array=self.line.split('\t')
        # scaffold5.1|size591573  AUGUSTUS    gene    1   21045   0.14    -   .   ID=symbB.v1.2.000001
        chr, item_type, start, stop, strand, attributes = line_array[0],line_array[2],line_array[3],line_array[4], line_array[6],line_array[8]
        if self.types_included:
            if not item_type in self.types_included:
                return
        for target_regexp in self.target_regexps:
            new_target_regexp = target_regexp + r'=([^;]+)'
            m = re.search(new_target_regexp,line_array[8])
            if not m.group(1):
                continue
            target = m.group(1)
            items_dict[target] = gff_items(start,stop,strand,chr)
        return items_dict


class gff_items():
    def __init__(self,start,stop,strand,chr):
        self.start=start
        self.stop=stop
        self.strand=strand
        self.chr=chr

