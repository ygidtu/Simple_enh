#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
2017.10.19

生成所有蛋白编码文件的bed类
"""
import os
from src.Bed import Bed
from src.Bed6 import Bed6

__dir__ = os.path.dirname(os.path.abspath(__file__))


class Gene:
    u"""
    返回符合要求的基因的区域
    """

    def __init__(self):
        u"""
        初始化，读取文件
        """
        beds = []
        with open(
                os.path.join(
                __dir__,
                'protein_coding.bed'
                )
        ) as reader:
            for line in reader:
                lines = line.split()
                if lines[-1] == '+':
                    start = int(lines[1]) - 1000
                    end = int(lines[2])
                else:
                    start = int(lines[1])
                    end = int(lines[2]) + 1000
                tem = Bed6(
                    chrom=lines[0],
                    start=start,
                    end=end
                )
                beds.append(tem)

        self.genes = Bed(beds)
