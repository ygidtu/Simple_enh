#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
后修改的Bed类，满足自己的需求.

2017.3.6 进一步进行修改和完善，补充了正负链的判断过程
2017.4.20 完成了对bed类的修改
2017.6.14 完成对脚本的拆分
共计两个Bed类

Bed6，储存每行Bed文件的类，支持bed3,6,12 以及非正式的4（位置信息，带正负链）
（1）包含正常bed6的各种setters和getters
        set_chrom (get_chrom): 设置（获取）染色体
        set_start (get_start): 设置（获取）起始位点
        set_end (get_end): 设置（获取）终止位点
        set_name (get_name): 设置（获取）名称信息
        set_score (get_score): 设置（获取）分值信息（该位置可以换成其他任意信息）
        set_strand (get_strand): 设置（获取）链信息
（2）重载了多个神奇方法
        __gt__、__eq__、__lt__、__ge__、__le__ 用来根据位置信息（染色体，起始终止位点，和正负链）判断Bed6的大小
        __hash__ 用来获取hash值，因此Bed6是iterable，支持in和not in（没有单独在重载__contain__等）等
        __and__ 用来判断Bed6之间是否具有重合位点
        __mul__ 用来判断Bed6的范围是否完全覆盖第二个Bed6
        __mod__ 用来计算两个Bed6的start site的距离
        __truediv__ 用来计算两个Bed6的中心位点的距离
（3）添加新方法
        cal_dist: 参数center（默认为False）, 用来根据参数，自动计算Bed6之间相应的距离
        __generate_chr_order: 自动生成chr1~chr22的染色体list，其他多余的chr会按照字符串的顺序导入
                              主要用于支持大小对比，并且用于Bed6的sort
        __get_index: 通过内置的__chrom来获取在染色体list的位置
"""
import sys


class Bed6():
    u"""Bed class."""

    __chrom = ''
    __start = 0
    __end = 0
    __name = ''
    __score = ''
    __strand = ''
    __append = ''  # 除了标准六行之外的信息
    __addition_chrom = []

    def __init__(self, chrom=None, start=None, end=None, name=None, score=None,
                 strand=None, append=None):
        u"""Init."""
        self.__order = self.__generate_chr_order()  # 必须先生成__order，随后需要用
        self.set_chrom(chrom)
        self.set_start(start)
        self.set_end(end)
        self.set_name(name)
        self.set_score(score)
        self.set_strand(strand)
        if append:
            if isinstance(append, list):
                append = '\t'.join(append)
            self.__append = append

    @staticmethod
    def __generate_chr_order():                     # 生成染色体的顺序
        order = []
        for i in range(1, 23):
            order.append('chr' + str(i))
        order.append('chr' + 'X')
        order.append('chr' + 'Y')
        order.append('chr' + 'M')
        return order

    # 各路setters
    def set_chrom(self, chrom):
        u"""Set chromosome."""
        if chrom not in self.__order:
            if chrom not in self.__addition_chrom:
                # 如果染色体不在默认名单里，就添加到额外的名单里
                self.__addition_chrom.append(chrom)
        self.__chrom = chrom

    def set_start(self, start):
        u"""Set start site."""
        if isinstance(start, int):
            self.__start = start
        else:
            self.__start = int(start)

    def set_end(self, end):
        u"""Set end site."""
        # 修改，完成，原始数据先转换成合适的类型，然后进行比对
        if isinstance(end, int):
            self.__end = end
        else:
            self.__end = int(end)
        if self.__end < self.__start:
            raise ValueError('end < start')

    def set_name(self, name):
        u"""Set track name."""
        if name is not None:
            self.__name = name

    def set_score(self, score):
        u"""Set track score, but also can set as other message."""
        if score is not None:   # 应该就是这出的问题
            self.__score = score

    def set_strand(self, strand):
        u"""Set track strand + or - or . ."""
        if strand is not None and strand in ('+', '-', '.'):
            self.__strand = strand

    def set_append(self, append):
        u"""Set append msg."""
        self.__append = append

    # 各种getters
    def get_chrom(self):
        u"""Return chromosome."""
        return self.__chrom

    def get_start(self):
        u"""Return track start site."""
        return self.__start

    def get_end(self):
        u"""Return track end site."""
        return self.__end

    def get_name(self):
        u"""Return track name."""
        return self.__name

    def get_score(self):
        u"""Return the message at column 5."""
        return self.__score

    def get_strand(self):
        u"""Return track strand."""
        return self.__strand

    def get_append(self):
        u"""Return append msg."""
        return self.__append

    def get_bed(self):
        u"""Return track."""
        import re
        standard = None
        for i in [self.__chrom, self.__start, self.__end, self.__name,
                  self.__score, self.__strand]:
            if i != '':                                     # 如果有一项为空，就不输出了
                if standard:
                    standard = standard + str(i) + '\t'
                else:
                    standard = str(i) + '\t'
        # 除去多余的\t，经测试，使用了正则才成功去除，rstrip和strip都失败，原因不明，可能我这不算字符串吧
        standard = re.sub('\t$', '', standard)
        # 原本standard和append之间还存在一个多余的\t，现在去除掉
        if self.__append:
            standard = standard + '\t' + self.__append
        return standard

    def get_center(self):
        u"""返回这个区域的中心."""
        return (self.__end - self.__start) / 2 + self.__start

    # 下面重载几个常用的选项，来满足自己的比较需求
    def __get_index(self, chrom):
        u"""这个get_index与Bed的不一样，是用来比较chrom谁比较靠前的."""
        try:
            if chrom in self.__order:              # 染色体正常就通过这个获取制定的排列顺序
                return self.__order.index(chrom)
            else:
                self.__addition_chrom.sort()       # 不常见染色体，就按字符串排序，意思意思
                self.__order.extend(self.__addition_chrom)
                return self.__order.index(chrom)
        except ValueError:
            print(chrom)
            sys.exit(1)

    def __lt__(self, other):                    # 大于, isupperstream
        u"""定义大于，用于比较两个bed track的先后循序，按照染色体，起始和终止位点， 上游."""
        if self.__get_index(self.__chrom) < self.__get_index(other.__chrom):
            return True
        elif self.__chrom == other.__chrom:
            if self.__end < other.__start:
                return True
        return False

    def __eq__(self, other):                # 相等
        u"""定义两条track相等的情况，这里只判断了链的位置信息，并不判断名称等等."""

        # 2017.8.18记，添加链的判断，和修复end位点的判断
        if self.__chrom == other.__chrom and \
                self.__start == other.__start and \
                self.__end == other.__end and \
                self.__compare_strand__(other):
            return True
        return False

    def __gt__(self, other):                # 小于  isdownstream
        u"""定义小于的情况, 下游."""
        if self.__get_index(self.__chrom) > self.__get_index(other.__chrom):
            return True
        elif self.__chrom == other.__chrom:
            if self.__start > other.__end:
                return True
        return False

    def __le__(self, other):
        u"""
        定义大于等于，用于比较两个bed track的先后循序，按照染色体，起始和终止位点.

        大于等于，两个区域可有可无交集，但是self的起始位点>=other，并且添加链的判断
        """
        if self.__get_index(self.__chrom) < self.__get_index(other.__chrom):
            return True
        elif self.__chrom == other.__chrom:
            if self.__start <= other.__start:
                if self.__compare_strand__(other) or \
                        (self.__strand == '+' or other.__strand == '-'):
                    return True
        return False

    def __ge__(self, other):            # 小于等于，两个区域可有可无交集，但是self的终止位点<=other
        u"""定义小于等于的情况."""
        if self.__get_index(self.__chrom) > self.__get_index(other.__chrom):
            return True
        elif self.__chrom == other.__chrom:
            if self.__end >= other.__end:
                if self.__compare_strand__(other) or \
                        (self.__strand == '-' or other.__strand == '+'):
                    return True
        return False

    def __compare_strand__(self, other):
        u"""用来判断链是否正确."""
        if self.__strand == other.__strand or\
                self.__strand in ('', '.') or\
                other.__strand in ('', '.'):
            return True
        return False

    def __and__(self, other):
        u"""用 & 来判断两个bed之间是否重合."""
        # 必须判断两个区域的染色体和正负链问题
        if self.__chrom == other.__chrom:
            if self.__compare_strand__(other):
                if self.__start <= other.__end and\
                        self.__end >= other.__start:
                    return True
                elif other.__start <= self.__end and\
                        other.__end >= self.__start:
                    return True
        return False

    def __hash__(self):
        u"""只有指定了hash方式，才能用于等于比较，用于in, not in等."""
        return hash(self.get_bed())

    def __mul__(self, other):
        u"""用来统计该bed是否包含有另外一个bed类的范围的."""
        if self.__chrom == other.__chrom and \
            self.__compare_strand__(other) and \
            self.__start <= other.__start and\
                self.__end >= other.__end:
            return True

        return False

    def __mod__(self, other):
        u"""取余 %，用来返回两个start site的距离."""
        return self.__start - other.__start

    def __truediv__(self, other):
        u"""除法/，用来返回两个bed中心的距离."""
        return self.get_center() - other.get_center()

    def isdownstream_same(self, other):
        u"""供downstream_target使用的判断上下游用的."""
        if self.__chrom == other.__chrom:
            if self.__start <= other.__end:
                return True
            else:
                return 2
        elif self.__get_index(self.__chrom) > self.__get_index(other.__chrom):
            return 2
        else:
            return 3

    def is_same_chrom(self, other):
        u"""
        用来判断两个bed6的染色体的大小
        """
        if self.__chrom == other.__chrom:
            return True
        elif self.__get_index(self.__chrom) > self.__get_index(other.__chrom):
            return 1
        else:
            return -1
