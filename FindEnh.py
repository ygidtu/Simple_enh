#!/usr/bin/python3
u"""
2017.3.8 用来通过自己的标准预测enhancer.

所需的是经过Process_Signal进行过归一化之后的结果

v1.0：最传统的三项甲基化修饰的enhancer预测和查找
V1.1：20170316 毕竟上次写完没测试，这次已经将get_hits调整好了，
      所以这次的主要目的就是，直接写好怎么提取相应的结果并且进行计算

v2.1: 重大失误，我的结果是没扩大上下2kb的，所以需要重新出一个扩大过的文件来进行比较文件来进行比对

V2.2: 添加了一个小函数，用来生成上游多少bp距离的蛋白编码文件 2017.4.23
      修改了使用标准化之后的bw文件，如果不是如此，比对就没有意义了
V3.0: 2017.5.17 周到enhancer后，添加了一个整合的过程

"""
import os
import re
import sys
import argparse
import textwrap
from shutil import rmtree
from datetime import datetime
from src.Bed import Bed
from src.Bed6 import Bed6
from FindEnh.bigFiles import BigBed
from FindEnh.bigFiles import BigWig
from FindEnh.Proteins import Gene

__version__ = "3.0"
__author__ = 'Yiming Zhang'
__email__ = "ygidtu@163.com"
__dir__ = os.path.dirname(os.path.abspath(__file__))


class FindEnhancer:
    u"""重新写，改成类，这里是几种enhancer组蛋白修饰的位点比较."""

    def __init__(self, indir=None, files=None):
        u"""初始化."""
        keys = ['Dnase-seq', 'H3K4me1', 'H3K4me3', 'H3K27ac']

        # 判断参数
        if indir is None:
            if isinstance(files, dict):
                if all([x for x in files.keys() if x in keys]):
                    self.files = files
                else:
                    sys.exit()
            else:
                sys.exit()
        else:
            self.files = self.__find_files__(indir)
        self.gene = Gene().genes
        self.tem_dir = os.path.join(__dir__, 'tem')

        if not os.path.exists(self.tem_dir):
            os.makedirs(self.tem_dir)
        pass

    @staticmethod
    def __find_files__(indir):
        u"""
        找到所有配对的bigbed和bigwig文件
        :param indir: 输入文件夹
        :return:
        """
        found = {}

        # 通过正则直接匹配出所需的细胞系，信号值和文件名
        base = r"^.*(\\|/)" \
               r"(?P<tissue>.*)(\\|/)" \
               r"(?P<signal>[\w-]+)(\\|/)" \
               r"(?P<file>\w+\.big(bed|wig))$"
        for parent, _, files in os.walk(indir):
            for file_ in files:
                path = os.path.join(parent, file_)

                match = re.search(base, path, re.I)
                if match:
                    match = match.groupdict()

                    # 复杂的构建字典的过程，此处可用递归改善
                    tissue = {
                        match['signal']: {
                            match['file'].split('.')[1]: path
                        }
                    }
                    if match['tissue'] in found.keys():
                        tissue = found[match['tissue']]

                        tem = {match['file'].split('.')[1]: path}
                        if match['signal'] in tissue.keys():
                            tem.update(tissue[match['signal']])
                        tissue.update({
                            match['signal']: tem
                        })

                    found.update({match['tissue']: tissue})
        return found

    def compare_peaks(self, files):
        u"""
        比对各种信号的区域值
        :return: 调整过区域范围的Bed
        """
        print("Start moving DHS region")

        dhs = files['Dnase-seq']['bigBed']
        dhs = BigBed(dhs)
        dhs.move_range(BigWig(files['Dnase-seq']['bigWig']))
        dhs = [x for x in dhs]

        print("Start comparing Peaks")
        signals = ['H3K27ac', 'H3K4me1']
        for signal in signals:
            if signal in files.keys():
                print("Comparing %s" % signal)
                bigbed = BigBed(files[signal]['bigBed'])
                bigwig = BigWig(files[signal]['bigWig'])
                peaks = bigwig.find_peaks(bigbed)
                dhs = self.__compare_peaks__(dhs, peaks)

        """
        最终在转化成Bed6类，方便后边处理
        这个地方设计复杂了，多做了一组转化，转过去在转回来
        """
        converted = []
        for bed in dhs:
            converted.append(Bed6(chrom=bed[0], start=bed[1], end=bed[2]))
        return converted

    @staticmethod
    def __compare_peaks__(base, peaks, gap=2000):
        i, j = 0, 0
        legal = set()
        base = sorted(base, key=lambda x: (x[0], x[1], x[2]))
        peaks = sorted(peaks, key=lambda x: (x[0], x[1]))
        while i < len(base) and j < len(peaks):
            bed = base[i]
            peak = peaks[j]
            if bed[0] < peak[0]:
                i += 1
            elif bed[0] > peak[0]:
                j += 1
            else:
                if peak[1] < bed[1] - gap:
                    j += 1
                elif peak[1] > bed[1] + gap:
                    i += 1
                else:
                    if bed[1] - gap <= peak[1] < bed[1] or\
                            bed[2] < peak[1] <= bed[2] + gap:
                        legal.add(bed)
                        i += 1
                    else:
                        j += 1
        return list(legal)

    def compare_signal(self, dhs, files, gap=2000):
        u"""
        符合要求的范围内的信号值
        :return:
        """
        print("Comparing Signal")
        legal = dhs
        if 'H3K4me1' in files.keys() and 'H3K4me3' in files.keys():
            legal = []

            me1 = BigWig(files['H3K4me1']['bigWig'])
            me3 = BigWig(files['H3K4me3']['bigWig'])

            for bed in dhs:
                if self.__compare_signal__(me1, bed, gap) >\
                        self.__compare_signal__(me3, bed, gap):
                    legal.append(bed)

        # 调整区域范围
        results = []
        for bed in legal:
            results.append(Bed6(
                chrom=bed.get_chrom(),
                start=bed.get_start() - gap,
                end=bed.get_end() + gap
            ))

        results = Bed(results) | self.gene
        return results

    @staticmethod
    def __compare_signal__(bigwig, peaks, gap=2000):
        u"""
        比对扩大过得范围内。不包括dhs范围本身内部的信号值
        :param bigwig:
        :param peaks:
        :return:
        """
        whole = bigwig.sum_peaks(
            peaks.get_chrom(),
            peaks.get_start() - gap,
            peaks.get_end() + gap
        )
        dhs = bigwig.sum_peaks(
            peaks.get_chrom(),
            peaks.get_start(),
            peaks.get_end()
        )
        return whole - dhs

    def find_enhancer(self):
        u"""
        找到符合要求的enhancer
        :return:
        """

        print("Start to find enhancer")
        enhancers = {}
        count = 1
        for tissue, files in self.files.items():
            print('#')
            begin = datetime.now()
            print('%d/%d %s' % (count, len(self.files), tissue))

            temfile = os.path.join(self.tem_dir, tissue)
            if os.path.exists(temfile):
                enhancer = Bed(temfile)
            else:
                peaks = self.compare_peaks(files)
                enhancer = Bed(self.compare_signal(peaks, files))
                enhancer.save(temfile)
            enhancers.update({tissue: enhancer})

            count += 1
            print("Time: %s" % str(datetime.now() - begin))
        rmtree(self.tem_dir)
        return enhancers


def main(outdir, indir=None, files=None):
    u"""
    开始找所有的enhancer
    :param indir: 输入文件夹
    :param outdir: 输出文件夹
    :param files: 输入的字典格式的文件路径
    :return: None
    """
    # 先找合适的enhancer
    if indir:
        enhancers = FindEnhancer(indir=indir).find_enhancer()
    elif files:
        enhancers = FindEnhancer(files=files).find_enhancer()
    else:
        sys.exit()

    for tissue, value in enhancers.items():
        value.save(os.path.join(outdir, tissue))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    subparser = parser.add_subparsers(
        help="batch run base on Directory architecture")
    group1 = subparser.add_parser(
        "batch",
        help="batch run",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
    ...         Please put files in right Directory architecture!
    ...         --------------------------------
    ...             Three level directories are needed!!
    ...             Top: cell or tissue names or just the label you want
    ...             Middle: Signal names, Must be one of Dnase-seq、H3K4me1、H3K4me3、H3K27ac
    ...             Bottom: a pair of bigBed and bigWig file
    ...         ''')
    )
    group1.add_argument(
        "-i",
        "--indir",
        help="Input directory, path/cells/signals/bigbed,bigwig"
    )
    group1.add_argument("-o", "--outdir", help="Output directory")

    group2 = subparser.add_parser(
        "single", help="single run, Just run for one set of files")
    group2.add_argument(
        "-l", "--label", help="default: enhancer", default="enhancer")
    group2.add_argument("-o", "--outdir", help="Output directory")
    group2.add_argument(
        "-dhsbb", help="path to Dnase-seq bigBed file", metavar="DHS bigBed")
    group2.add_argument(
        "-dhsbw", help="path to Dnase-seq bigWig file", metavar="DHS bigWig")
    group2.add_argument(
        "-27acbb", help="path to H3K27ac bigBed file", metavar="H3K27ac bigBed")
    group2.add_argument(
        "-27acbw", help="path to H3K27ac bigWig file", metavar="H3K27ac bigWig")
    group2.add_argument(
        "-me1bb", help="optional, path to H3K4me1 bigBed file", metavar="H3K4me1 bigBed")
    group2.add_argument(
        "-me1bw", help="optional, path to H3K4me1 bigWig file", metavar="H3K4me1 bigWig")
    group2.add_argument(
        "-me3bw", help="optional, path to H3K4me3 bigWig file", metavar="H3K4me3 bigWig")

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = vars(parser.parse_args())
        if 'indir' in args.keys():
            # 如果没有指定任何参数，输出帮助信息
            if any([True for x in args.values() if x is None]):
                parser.print_help()
                sys.exit()
            main(indir=args['indir'], outdir=args['outdir'])
        elif 'dhsbb' in args.keys():
            files = {}
            label = args.pop('label')
            outdir = args.pop('outdir')

            if outdir is None:
                print("Please set outdir")
                parser.print_help()
                sys.exit()

            __para__ = {
                'dhs': 'Dnase-seq',
                '27ac': 'H3K27ac',
                'me1': 'H3K4me1',
                'me3': 'H3K4me3'
            }

            # 将指定的文件整理成需要的结果
            for key, value in args.items():
                if value is None:
                    continue
                if key.endswith('bb'):
                    tem = {'bigBed': value}
                elif key.endswith('bw'):
                    tem = {'bigWig': value}

                key = __para__[key[:-2]]
                if key in files.keys():
                    tem.update(files[key])
                files.update({key: tem})

            if len(files) < 1:
                print("Please set bigBed and bigWig files")
                parser.print_help()
                sys.exit()

            files = {label: files}
            main(outdir=outdir, files=files)
