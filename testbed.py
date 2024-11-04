import thermo

fug = thermo.eos_mix.PRMIX([190.564, 305.322], [4.5992E6, 4.8722E6], [0.01142, 0.0995], [0.5, 0.5], [[0, -0.003], [-0.003, 0]], T=150, P=.05*1E6)

print(fug.fugacities_g)