# Interactions between bodyparts

When parts of different kinds collide with each other various things happen:

### CELL:
- CELL - nothing
- SPIKE - cell gets damaged
- SHIELD - nothing
- BOUNCE - cell gets pushed back

### SPIKE:
- CELL - cell gets damaged
- SPIKE - both spikes get damaged. 
  But also in the current implementation this can't happen, they'd cross
  each other and collide with the cells they are attached to.
- SHIELD - nothing
- BOUNCE - bounce is deflated, spike gets pushed back strongly

### SHIELD:
- CELL - nothing
- SPIKE - nothing
- SHIELD - nothing
- BOUNCE - shield gets pushed back

### BOUNCE:
- CELL - cell gets pushed back
- SPIKE - bounce is deflated, spike gets pushed back strongly
- SHIELD - shield gets pushed back
- BOUNCE - undecided.
  - bounce off normally
  - bounce off strongly but do not deflate
  - nothing, bounce cancel each other out
  - something weird, like transferring momentum so both users travel together
