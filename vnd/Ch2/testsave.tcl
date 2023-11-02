# Visual Neuronal Dynamics
# Visualization State
::NeuronVND::loadFiles /Scr/mariano/Tutorials/VND/Ch2/circuit_config.json false false
# List of representations
::NeuronVND::createRepArgs show true style soma color red material Opaque selection {soma x < 0}
::NeuronVND::createRepArgs show true style morphology_spheretube color Type material Opaque selection {stride 100}
::NeuronVND::createRepArgs show true style soma color Type material Opaque selection {(soma x > 0) && (soma y > 0)}
