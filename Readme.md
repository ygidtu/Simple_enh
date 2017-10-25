Find enhancers
===
## pyBigWig are required, `pip3 install pyBigWig`
pybigWig: [<https://github.com/dpryan79/pyBigWig>](https://github.com/dpryan79/pyBigWig)


**Base on Dnase-seq data and 3 histone modification data**
- H3K27ac
- H3K4me1
- H3K4me3

bigBed format peaks files & bigWig format signal files are required

#### Two diffent ways to use this script


#### 1. single
> the bigBed and bigWig files of one signal should paired.
```bash
python3 Find_enhancer.py single -h

usage: Find_enhancer.py single [-h] [-l LABEL] [-o OUTDIR] [-dhsbb DHS bigBed]
                               [-dhsbw DHS bigWig] [-27acbb H3K27ac bigBed]
                               [-27acbw H3K27ac bigWig]
                               [-me1bb H3K4me1 bigBed] [-me1bw H3K4me1 bigWig]
                               [-me3bw H3K4me3 bigWig]

optional arguments:
  -h, --help            show this help message and exit
  -l LABEL, --label LABEL
                        default: enhancer
  -o OUTDIR, --outdir OUTDIR
                        Output directory
  -dhsbb DHS bigBed     path to Dnase-seq bigBed file
  -dhsbw DHS bigWig     path to Dnase-seq bigWig file
  -27acbb H3K27ac bigBed
                        path to H3K27ac bigBed file
  -27acbw H3K27ac bigWig
                        path to H3K27ac bigWig file
  -me1bb H3K4me1 bigBed
                        optional, path to H3K4me1 bigBed file
  -me1bw H3K4me1 bigWig
                        optional, path to H3K4me1 bigWig file
  -me3bw H3K4me3 bigWig
                        optional, path to H3K4me3 bigWig file
```

### batch
```bash
python3 Find_enhancer.py batch -h

usage: Find_enhancer.py batch [-h] [-i INDIR] [-o OUTDIR]

...         Please put files in right Directory architecture!
...         --------------------------------
...             Three level directories are needed!!
...             Top: cell or tissue names or just the label you want
...             Middle: Signal names, Must be one of Dnase-seq、H3K4me1、H3K4me3、H3K27ac
...             Bottom: a pair of bigBed and bigWig file
...         

optional arguments:
  -h, --help            show this help message and exit
  -i INDIR, --indir INDIR
                        Input directory, path/cells/signals/bigbed,bigwig
  -o OUTDIR, --outdir OUTDIR
                        Output directory
```

for example, the data files of HepG2 should be placed like below:

- xxx/HepG2/Dnase-seq/
    - xxx.bigBed
    - xxx.bigWig
- xxx/HepG2/H3K27ac/
    - xxx.bigBed
    - xxx.bigWig

then you could just run: 
- `pyhton3 Find_enhancer.py batch -i xxx/ -o test`
- `python3 Find_enhancer.py batch -i xxx/HepG2/ -o test` 
    - just calculate the enhancer of HepG2

#### Please note that do not please same format of files in same direcotry (2 bigBed files in same directory), and do not name those files with special symbol

