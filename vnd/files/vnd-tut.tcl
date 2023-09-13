# Visual Neuronal Dynamics
# Visualization State
::NeuronVND::loadFiles /home/student/Tutorials/vnd/files/300_cells/circuit_config.json false false
# List of representations
::NeuronVND::createRepArgs show true style soma color silver material Opaque selection {soma x < 0}
::NeuronVND::createRepArgs show true style soma color red material Opaque selection {(soma x < 0) && (soma y > 0)}
::NeuronVND::createRepArgs show true style morphology color yellow material Opaque selection {stride 100}
