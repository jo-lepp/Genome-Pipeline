
process name_here {

container 'container_name' //using docker? how do that? -nextflow.config?

input:
val variable
path input_files

output:
path '${variable}-output', emit: short_name // i think this is what the emit tag is for?
path 'path' , emit: second_out

script:
"""
// put commands here
"""

}

include {process_name} from 'path_to_module'

params{ //move to nextflow.config, have people enter in their values before running? (param file)
// default values for inputs?
// i'm confused about this bit
input: String = "default value"
path_input: String = "/home/jalepper/Genome-Pipeline"
}

workflow {

main:
channel_name = channel.of(can be a bunch of inputs)
name_here(params.input or channel_name)
    .map() 

// using the output of one process in another process
second_process(name_here.out)

channel_two = channel.fromPath(params.path_input)

}

output {

// what's the point of this?
first_output {
    path = 'path here'
    // are paths absolute?
    mode = 'copy'
}

}

//params file would be so helpful