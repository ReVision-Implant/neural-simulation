# Visual Neuronal Dynamics
# Visualization State
# List of representations
::NeuronVND::createRepArgs show false style soma color red material Opaque selection {soma x < 0}
::NeuronVND::createRepArgs show true style morphology color Type material Opaque selection {stride 100}
::NeuronVND::createRepArgs show true style soma color Type material Opaque selection {(soma x > 0) && (soma y > 0)}
