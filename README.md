# Genome-Pipeline
Work in progress, creating an apptainer script to run genome assembly

# Cloning the GitHub Repository

```bash
git clone https://github.com/jo-lepp/Genome-Pipeline.git
```
# Apptainer Container setup

This setup only needs to be done once. Navigate to the cloned github folder and run the following commands.

```bash
apptainer build env.sif env.def
apptainer build env_fast.sif env_fast.def

apptainer build fast_mode.sif fast_mode.def
apptainer build full_mode.sif full_mode.def
apptainer build pre_assembled.sif pre_assembled.def
```

# Full Mode

Runs: Myloasm, MetaBat2, Remag, Semibin (soil, self, global), dRep, Busco, Checkm2, gtdbtk, gvclass

Command: 
```bash
apptainer run --bind /home/user/path_to_data_folder:/home/user/path_to_data_folder --bind /home/user/path_to_output_folder:/home/user/path_to_output_folder --env input_fgz=/home/user/path/to/fgz --env output_path=/home/user/path/to/output_folder /path/to/folder/full_mode.sif
```

Explanation:
–bind: Connect specific folders so that the container can view what’s inside and use those files. Path_to_data_folder should be the folder that contains the fastq.gz, path_to_output_folder should be the folder where you want the outputs to land.

–env: setting variables for the container.

/path/to/folder should be the file path to the cloned github folder.

# Fast Mode

Runs: Myloasm, Metabat2, Semibin (soil/global), dRep, busco, checkm2, cmsearch (euk)
Doesn't Run: Remag, Semibin self-trained model, gvclass, gtdbtk

Command:
```bash
apptainer run --bind /home/user/path_to_data_folder:/home/user/path_to_data_folder --bind /home/user/path_to_output_folder:/home/user/path_to_output_folder --env input_fgz=/home/user/path/to/fgz --env output_path=/home/user/path/to/output_folder /path/to/folder/fast_mode.sif
```

Explanation:
–bind: Connect specific folders so that the container can view what’s inside and use those files. Path_to_data_folder should be the folder that contains the fastq.gz, path_to_output_folder should be the folder where you want the outputs to land.

–env: setting variables for the container.

/path/to/folder should be the file path to the cloned github folder.

# Pre-assembled

Runs: MetaBat2, Remag, Semibin (soil, self, global), dRep, Busco, Checkm2, gtdbtk, gvclass

Essentially, the fast-mode version with a given assembly file

Command:
```bash
apptainer run --bind /home/user/path_to_data_folder:/home/user/path_to_data_folder --bind /home/user/path_to_output_folder:/home/user/path_to_output_folder --bind /home/user/path_to_assembly_folder:/home/user/path_to_assembly_folder: --env input_fgz=/home/user/path/to/fgz --env output_path=/home/user/path/to/output_folder --env assembly=/home/user/path/to/assembly.fa /path/tofolder/pre_assembled.sif
```

Explanation:
–bind: Connect specific folders so that the container can view what’s inside and use those files. Path_to_data_folder should be the folder that contains the fastq.gz, path_to_output_folder should be the folder where you want the outputs to land. path_to_assembly_folder should contain the assembly. Note: if the assembly and fastq.gz are in the same folder, you only need one bind statement.

–env: setting variables for the container. 

/path/to/folder should be the file path to the cloned github folder.




