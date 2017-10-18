#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
2017.10.15
通过三种信号的原始文件来找相应的enhancer
"""
import os
import sys
import argparse
from subprocess import check_call
from subprocess import PIPE
from subprocess import Popen
from subprocess import CalledProcessError
from subprocess import STDOUT
try:
    import pyBigWig
except ImportError as err:
    print("Please install pyBigWig first")
    print("pip3 install pyBigWig")
    sys.exit()

try:
    null = open(os.devnull, 'w')
    check_call('bwtool --version', shell=True, stdout=null, stderr=STDOUT)
    null.close()
except CalledProcessError:
    print('Please install bwtool properly first')
    sys.exit()
from Bed import Bed
from Bed6 import Bed6


__dir__ = os.path.dirname(os.path.abspath(__file__))


class DealPeaks:
    u"""
    处理甲基化信号
    """

    def __init__(self, args):
        u"""
        初始化，提供文件目录
        """
        self.args = vars(args)
        self._check_file_()
        self.dhs = self._compare_()

    def _check_file_(self):
        u"""
        检查文件合法性
        """
        for key, value in self.args.items():
            if key.endswith('bb'):
                if value is not None and os.path.exists(value):
                    bb = pyBigWig.open(value)
                    if bb.isBigBed():
                        bb.close()
                        continue
                raise ValueError('%s - %s is not a bigBed file' % (key, value))

    def _readbb_(self, infile):
        u"""
        读取bb文件
        """
        print('Reading:\t%s' % infile)
        bb = pyBigWig.open(infile)

        results = []
        for chrom, size in bb.chroms().items():
            try:
                data = bb.entries(chrom, 0, size)
            except RuntimeError:
                print(chrom)
                continue
            for tem in data:
                new_bed = Bed6(
                    chrom=chrom,
                    start=tem[0],
                    end=tem[1]
                )
                results.append(new_bed)

        return Bed(results)

    def _compare_(self):
        u"""
        比对dhs区域
        """
        dhs = self._readbb_(self.args['dhsbb'])
        me1 = self._readbb_(self.args['me1bb'])
        ac27 = self._readbb_(self.args['27bb'])

        print('Comparing DHS & H3K4me1')
        tem = dhs & me1
        tem = Bed(list(tem.keys()))
        print('Comparing DHS & H3K27ac')
        tem & ac27
        return tem

    def save(self, outfile):
        u"""
        保存文件
        """
        print('Comparing results are saving to:\t%s' % outfile)

        self.dhs.save(outfile)


class DealSignal:
    u"""
    比对每个区域内的信号
    """

    def __init__(self, bed, args):
        u"""
        初始化
        """
        self.bed = bed
        self.me1, self.me3 = self._check_bw_files_(vars(args))

    @staticmethod
    def _check_bw_files_(args):
        u"""
        检查文件合法性
        """
        files = []
        for i in ['me1bg', 'me3bg']:
            path = args.get(i)
            if path is not None and os.path.exists(path):
                bw = pyBigWig.open(path)
                if bw.isBigWig:
                    files.append(path)
                else:
                    raise ValueError('%s - %s is not BigWig file' % (i, path))
            else:
                raise ValueError('%s - %s is not BigWig file' % (i, path))
        return tuple(files)

    @staticmethod
    def _call_bwtool_(bed, bigwig):
        u"""
        调用bwtool
        """
        print('bwtool start processing %s' % bigwig)
        commands = [
            'bwtool',
            'summary',
            bed,
            bigwig,
            'stdout',
            '-with-sum'
        ]
        data = {}
        try:
            proc = Popen(commands, stdout=PIPE)

            while True:
                line = proc.stdout.readline().decode('utf-8')

                if not line:
                    break

                lines = line.split()

                key = Bed6(
                    chrom=lines[0], start=int(
                        lines[1]), end=int(lines[2])
                )
                data.update({key: float(lines[-1])})
        except CalledProcessError as err:
            print(err)
            sys.exit()
        return data

    def save(self, outfile):
        u"""
        保存结果
        """
        me1 = self._call_bwtool_(self.bed, self.me1)
        me3 = self._call_bwtool_(self.bed, self.me3)

        print('Start final comparing')
        with open(outfile, 'w+') as writer:
            for key in me1.keys():
                if me1[key] > me3[key]:
                    writer.write(key.get_bed() + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="利用DHS高敏位点及三种信号文件寻找enhancer")
    parser.add_argument(
        '-dhsbb',
        help="DHS的bigBed文件, bigBed file of Dnase-seq"
    )
    parser.add_argument(
        '-me1bb',
        help="H3K4me1的bigBed文件, bigBed file of H3K4me1"
    )
    parser.add_argument(
        '-27bb',
        help="H3K27ac的bigBed文件, bigBed file of H3K27ac"
    )
    parser.add_argument(
        '-me1bg',
        help="H3K4me1的bigWig文件, bigWig file of H3K4me1"
    )
    parser.add_argument(
        '-me3bg',
        help="H3K4me3的bigWig文件， bigwig file of H3K4me3"
    )

    parser.add_argument(
        '-outfile',
        help="保存最终结果，results files",
        default='enhancer.bed'
    )

    if len(sys.argv) == 1:
        parser.print_help()
    outfile = parser.parse_args().outfile
    DealPeaks(parser.parse_args()).save(outfile + '.tem')
    DealSignal(outfile + '.tem', parser.parse_args()).save(outfile)
    os.remove(outfile + '.tem')
    print('All processes are finished')
