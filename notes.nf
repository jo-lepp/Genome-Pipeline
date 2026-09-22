
process name_here {

input:
val variable

output:
path '${variable}-output'

script:
"""
// put commands here
"""

}

params{
// default values for inputs?
input: String = "default value"
path_input: String = "/home/jalepper/Genome-Pipeline"
}

workflow {

main:
channel_name = channel.of(can be a bunch of inputs)
name_here(params.input or channel_name)
    .map() // might be good if we decide to pull best bins

channel_two = channel.fromPath(params.path_input)


publish:

// outputs here

}

output {

// what's the point of this?

}