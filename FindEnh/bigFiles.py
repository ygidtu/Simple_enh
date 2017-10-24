#!/usr/bin/env python3
u"""
2017.10.19

重新订立标准之后，重新进行分析的
定义bigbed类
"""
import os
import sys
try:
    import pyBigWig
except ImportError as err:
    print(err)
    sys.exit(err)
from src.Bed6 import Bed6


class BaseBig:
    u"""
    bigBed和bigWig的基础类
    """

    def __init__(self, infile):
        u"""
        初始化
        """
        self.bigbed, self.bigwig = None, None
        if os.path.exists(infile):
            bigfile = pyBigWig.open(infile)

            if bigfile.isBigWig():
                self.bigwig = bigfile
            elif bigfile.isBigBed():
                self.bigbed = bigfile
            else:
                raise ValueError('%s is not a legit file' % infile)
            self.chroms = bigfile.chroms()

    @staticmethod
    def __find_where_max__(intervals):
        u"""
        返回某段信号内，峰值最高的点
        :param intervals: 
        :return: 
        """
        index = 0   # 记录index
        where = 0   # 记录最高点
        now = None  # 记录最大值
        for i in intervals:
            if now is None:
                now = i
            elif i > now:
                now = i
                where = index

            index += 1
        return where


class BigBed(BaseBig):
    u"""
    读取bigbed的类，并且做操作的类
    """

    def __init__(self, infile):
        u"""
        初始化
        """
        BaseBig.__init__(self, infile)
        if self.bigbed is None:
            raise ValueError('%s is not bigbed file' % infile)

        self.bed = self.__read__()

    def __iter__(self):
        u"""
        重载iter类
        :return: bigbed文件中的原始值
        """
        for bed in self.bed:
            yield bed.get_chrom(), bed.get_start(), bed.get_end()

    def __read__(self):
        u"""
        从bigbed文件中，读取所有的位点信息
        """
        bed = []
        for chrom, length in self.chroms.items():
            for data in self.bigbed.entries(chrom, 0, length):
                # pyBigWig返回的值中不包含chrom的信息，此处手动添加
                tem = Bed6(
                    chrom=chrom,
                    start=data[0],
                    end=data[1]
                )
                bed.append(tem)

        return sorted(bed)

    def move_range(self, bigwig):
        u"""
        迁移中心区域
        :param bigwig: 
        :return: 
        """
        if not isinstance(bigwig, BigWig):
            raise ValueError('Please submit BigWig class')

        chroms = bigwig.bigwig.chroms()
        not_match = set()
        new_bed = []
        for bed in self.bed:

            if bed.get_chrom() not in chroms.keys() or\
                    bed.get_end() > chroms[bed.get_chrom()]:
                not_match.add(bed.get_chrom())
                continue

            values = bigwig.bigwig.values(
                bed.get_chrom(),
                bed.get_start(),
                bed.get_end()
            )
            # 计算中心点
            center = self.__find_where_max__(values) + bed.get_start()
            # 计算中心周边的距离
            length = (bed.get_end() - bed.get_start()) // 2
            # 计算新的迁移过的区域
            bed.set_start(center - length)
            bed.set_end(center + length)

            # 返回新的
            new_bed.append(bed)
        if len(not_match) > 0:
            print("bigBed and bigWig files are not matching")
            print(sorted(list(not_match)))
        self.bed = new_bed

    def __len__(self):
        u"""
        重载长度
        :return: 
        """
        return len(self.bed)

    def beds(self):
        u"""
        generator for beds
        :return: 
        """
        for bed in self.bed:
            yield bed


class BigWig(BaseBig):
    u"""
    处理bigwig类
    """

    def __init__(self, infile):
        u"""
        初始化
        :param infile: bigwig files
        """
        BaseBig.__init__(self, infile)
        if self.bigwig is None:
            raise ValueError('%s is not bigwig file' % infile)

    def find_peaks(self, bigbed):
        u"""
        记录最高点的位置，与信号值
        :param bigbed: 对应的bigwig的peaks位点
        :return: 
        """
        if not isinstance(bigbed, BigBed):
            raise ValueError('Please submit BigWig class')

        peaks = []
        chroms = self.bigwig.chroms().keys()
        not_match = set()
        for bed in bigbed:
            if bed[0] not in chroms:
                not_match.add(bed[0])
                continue
            values = self.bigwig.values(
                bed[0],
                bed[1],
                bed[2]
            )
            where_max = self.__find_where_max__(values)

            tem = [bed[0], bed[1] + where_max]

            peaks.append(tem)

        if len(not_match) > 0:
            print("bigBed and bigWig files are not matching")
            print(sorted(list(not_match)))
        return peaks

    def sum_peaks(self, *args):
        u"""
        统计某区域内的信号值
        :param peaks: [chr, start, end]
        :return: sum
        """

        values = self.bigwig.values(
            args[0],
            args[1],
            args[2]
        )
        return sum(values)


if __name__ == '__main__':
    pass
