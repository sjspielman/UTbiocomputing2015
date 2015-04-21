# TACC Basics

Information about TACC for bioinformatics is sparse on the web; the BioITeam Wiki (search for it) has some info. TACC has user guides, and a *very* useful consulting system for users. You can also email me or anyone else at the Bioinformatics Consulting Group.

### TACC Website
`portal.tacc.utexas.edu`

- User Guides
- Consulting Ticket System. I've been impressed with their speed in answering questions (even very basic ones!), and their willingness to help me with odd problems.


### Disk Sizes
- $HOME (`cdh` shortcut)
	- 5GB on Stampede, 1GB on Lonestar
	- Backed-up
- $WORK (`cdw` shortcut)
	- 400GB on Stampede, 250GB on Lonestar
	- Not backed-up
- $SCRATCH (`cds` shortcut)
	- 8.5PB (8,500,000GB) on Stampede, 1.4PB on Lonestar
	- Not backed-up, all files older than 10 days are deleted

### Modules
- Search for a module: `module spider <search_term>`
- List current modules: `module list`
- Load a module: `module load <module_name>`
- Swap modules: `module swap <old_module> <new_module>`

### Queues
- `normal`: best for most job submissions.
	- Stampede: 48 hour limit, 16 cores per node
	- Lonestar: 24 hour limit, 12 cores per node
- `largemem`: when you need *lots* of memory, these have 1TB.
	- Stampede: 48 hour limit, 32 cores per node
	- Lonestar: 24 hour limit, 24 cores per node
- `development`: good for quick, small tests.
	- Stampede: 2 hour limit, 16 cores per node
	- Lonestar: 1 hour limit, 12 cores per node

### Launcher Scripts
I recommend using `launcher_creator.py`. Information is available here:
https://wikis.utexas.edu/display/bioiteam/launcher_creator.py
(or you can search for "BioITeam wiki" and look under "Software")

To use `launcher_creator.py`, you should probably add a line to your `.profile` file (see last section).

Lonestar uses SGE as its job scheduling system, Stampede uses SLURM.

- Submit a job:
	- Stampede: `sbatch job.slurm`
	- Lonestar: `qsub job.sge`
- Check status:
	- Stampede: `qstat`
	- Lonestar: `qstat`
- Delete a job:
	- Stampede: `scancel job_id`
	- Lonestar: `qdel job_id`
	- (get `job_id` from `qstat`)

### Pre-rolled BioITeam Scripts
- BWA-MEM: `map_BWAmem`
- BLAST+: `split_blast`
- Trinity: `assemble_trinity`

### Add BioITeam To Your $PATH
One way to do it (but keep it on one line):
		
	echo 'source /corral-repl/utexas/BioITeam/bin/profile_ngs_course.bash' >> ~/.profile_user

# More Nodes Doesn't Make It Faster
