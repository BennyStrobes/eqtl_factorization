"""
MIT License
Copyright (c) 2020 Warren W. Kretzschmar
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from pathlib import Path

import pysam
import pdb


class BamWriter:
    def __init__(self, alignment, barcodes, prefix):
        self.alignment = alignment
        self.prefix = prefix
        self.barcodes = set(barcodes)
        self._out_files = {}

    def write_record_to_barcode(self, rec, barcode):
        if barcode not in self.barcodes:
            return
        if barcode not in self._out_files:
            self._open_file_for_barcode(barcode)
        self._out_files[barcode].write(rec)

    def _open_file_for_barcode(self, barcode):
        self._out_files[barcode] = pysam.AlignmentFile(
            f"{self.prefix}_{barcode}.bam", "wb", template=self.alignment
        )


def main(input_bam, barcodes_file, output_prefix, contigs):
    """Split a 10x barcoded sequencing file into barcode-specific BAMs

    input:
    barcodes_file: a file containing barcodes or a single barcode
    contigs: '.' for all contigs, 'chr1' for the contig 'chr1',
    or '1-5' for chromosomes 1, 2, 3, 4, and 5
    """
    alignment = pysam.AlignmentFile(input_bam)
    if Path(barcodes_file).is_file():
        with open(barcodes_file, "r") as fh:
            barcodes = [l.rstrip() for l in fh.readlines()]
    else:
        barcodes = [barcodes_file]
        print(f"Extracting single barcode: {barcodes}")
    writer = BamWriter(alignment=alignment, barcodes=barcodes, prefix=output_prefix)
    if contigs == ".":
        print("Extracting reads from all contigs")
        recs = alignment.fetch()
    else:
        if "-" in contigs:
            start, end = contigs.split("-")
            print(f"Extracting reads from contigs {start} to {end}")
            recs = (alignment.fetch(str(contig)) for contig in range(start, end + 1))
        elif "," in contigs:
            contigs = contigs.split(",")
            print(f"Extracting reads from contigs {contigs}")
            recs = (alignment.fetch(str(contig)) for contig in contigs)
        else:
            print("Extracting reads for one contig: {contigs}")
            recs = (alignment.fetch(c) for c in [contigs])
    for rec in recs:
        try:
            barcode = rec.get_tag("CB")
            writer.write_record_to_barcode(rec=rec, barcode=barcode)
        except KeyError:
            pass


if __name__ == "__main__":
    import sys

    main(
        input_bam=sys.argv[1],
        barcodes_file=sys.argv[2],
        output_prefix=sys.argv[3],
        contigs=sys.argv[4],
    )
