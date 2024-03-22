clear; cd src; make; cd ..
# src/slugs --swarmTest examples/DARS_Example_1.slugsin  DARS_Example_1_output
# src/slugs --swarmTest examples/simple_swarm_interactive_v2.slugsin  simple_swarm_interactive_v2_output
src/slugs --swarmTest examples/DARS_Demo_V1.slugsin  dars_demo_v1_output
. makePngs.sh
