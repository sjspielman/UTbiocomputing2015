# Introduction to TACC
Learn what a cluster is, and how to use the world-class clusters available at the Texas Advanced Computing Center. The course will discuss the basic architecture of the Lonestar and Stampede computing clusters, how they compare to a regular computer, job launchers and job scheduling, and how to submit your own jobs to TACC. Custom tools by the Bioinformatics Consulting Group for job submission will be emphasized. Basic familiarity with the UNIX command line will be assumed.

## Intro To TACC Slideshow
Quick summary of systems.

## "Supercomputer" vs "cluster" vs computer
Stampede is not like an ordinary powerful single computer, just 6,000 times faster. The fact that Stampede consists of 6,000 separate computers makes "supercomputer" a bit of a misnomer—"cluster" or "supercluster" might be more accurate. (Doesn't really matter though.) Large clusters are inherently "parallel". To get the benefit of running 1,000s of computers at once, you need to keep them all fed, and to keep them all fed, things need to run in parallel.

Let's compare the basic architectures of a single computer vs a cluster to get a better idea of what's going on here.


## Disk

Lots of users on Stampede, Lonestar. Each user needs space for permanent files (like scripts and code), large files (permanent data), and temporary files for computation.

- $HOME
	- 5GB on Stampede, 1GB on Lonestar
	- Backed-up, no purge
	- `cdh` shortcut
- $WORK
	- 400GB on Stampede, 250GB on Lonestar
	- Not backed-up
	- `cdw` shortcut
- $SCRATCH
	- 8.5PB (8,500,000GB) on Stampede, 1.4PB on Lonestar
	- Not backed-up, all files older than 10 days are deleted
	- TACC may delete files sooner than 10 days if they feel like it, or you looked at them funny
	- `cds` shortcut

## Queues

The overwhelming number of nodes on Lonestar and Stampede are identical. But there are some special-purpose nodes. The different types of nodes are organzied into "queues". The most common queue is the "normal" queue. (Show the info in the TACC User Guide.)

- Normal
	- 12 cores (Lonestar), 16 cores (Stampede)
	- 24GB RAM (Lonestar), 32GB RAM (Stampede)
	- 24 hrs (Lonestar), 48 hrs (Stampede)
	- 1 SU per core per hour
- Largemem
	- 24 cores (Lonestar), 32 cores (Stampede)
	- 1TB RAM
	- 24 hrs (Lonestar), 48 hrs (Stampede)
	- 4 SU (Lonestar), 2 SU (Stampede)
- Development
	- 1 hr (Lonestar), 2 hrs (Stampede)

## First Launcher
We'll create a launcher script on Lonestar for a simple task (just to print "Hello World"), and then take a look inside the launcher script. Scott, myself and John Hawkins have written a small program to make creating launchers a little easier.

- `launcher_creator.py -h`
- `launcher_creator.py -n hello_world -t 00:01:00 -a DNAdenovo -q normal -w 1 -b 'echo "hello world"'`
- `less hello_world.sge`
- `qsub hello_world.sge`
- `qstat`

Look at the .o and .e files, explain how they're connected to STDOUT and STDERR I mentioned last class (intro to Unix).

`qsub` is different on Stampede (you use `sbatch`), but `qstat` works.

## Modules
Lots of users in diverse want to use lots of programs with different dependencies and libraries. Different programs require different, often conflicting, compute environments. TACC solution is modules. (This is TACC-specific.)

- `module load bedtools`
- `module swap intel gcc`
- `module load trinityrnaseq`
- `module load bowtie`
- `module spider trinityrnaseq`

## Allocations
The TACC clusters are Free for UT-affiliated people. Other people have to wait in line.

- SU funny money.
- Allocation determines group
- 1-1 correspondence between allocations and TACC Unix groups
- Show TACC portal for allocations
- Show CCBB-TACC-Fall2014, Scott PI, me delegate
- Show consulting ticket system, user guides

## Exploiting Parallelism
The previous launcher that just printed "Hello world" was a standard starting point, but it didn't really take advantage of Lonestar's computational power. The program ran for less than a second on a single machine. Using a cluster effectively means using many nodes at once. But assigning more computers to a program doesn't necessarily make it faster. The program needs to be able to spread the work between multiple nodes.

## Paper Notes On Parallelism (Part B)

## Parallel BWA Example

- `mkdir BWA_Example`
- `cd BWA_Example`
- `cp /work/01863/benni/IntroToTacc/SRR* .`
- `split -d -l 1000 SRR1580546_1.fastq reads_R1-`
- `split -d -l 1000 SRR1580546_2.fastq reads_R2-`
- `nano bwa_commands.txt`

		bwa mem /corral-repl/utexas/BioITeam/tmp/benni/hg19_ref_index/Homo_sapiens.GRCh37.60.dna.fasta reads_R1-00 reads_R2-00 > test00.sam
		bwa mem /corral-repl/utexas/BioITeam/tmp/benni/hg19_ref_index/Homo_sapiens.GRCh37.60.dna.fasta reads_R1-01 reads_R2-01 > test01.sam
		bwa mem /corral-repl/utexas/BioITeam/tmp/benni/hg19_ref_index/Homo_sapiens.GRCh37.60.dna.fasta reads_R1-02 reads_R2-02 > test02.sam
		bwa mem /corral-repl/utexas/BioITeam/tmp/benni/hg19_ref_index/Homo_sapiens.GRCh37.60.dna.fasta reads_R1-03 reads_R2-03 > test03.sam

- `launcher_creator.py -a DNAdenovo -n bwa_test -t 00:10:00 -q normal -m 'module load bwa/0.7.7' -w 2 -j bwa_commands.txt`
- `qsub bwa_test.sge`
- `qstat`
- `awk '! /^@/' test02.sam > test_noheader-02.sam`

- Splitting files
	- `ubersplit.py`
	- Unix `split` command

## Custom BCG Scripts
- `map_BWAmem`
- `split_blast`
- `assemble_trinity`

## BioITeam $PATH
`source /corral-repl/utexas/BioITeam/bin/profile_ngs_course.bash` into `~/.profile_user`