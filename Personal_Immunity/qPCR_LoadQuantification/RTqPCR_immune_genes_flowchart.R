pacman::p_load(
  DiagrammeR, # for flow diagrams
  networkD3, # For alluvial/Sankey diagrams
  tidyverse
) # data management and visualization




DiagrammeR::grViz("               # All instructions are within a large character string
digraph surveillance_diagram {    # 'digraph' means 'directional graph', then the graph name 
  
  # graph statement
  #################
  graph [layout = dot,
         rankdir = TB,            # layout top-to-bottom
         fontsize = 10]
  

  # nodes (circles)
  #################
  node [shape = box,           # shape = circle
       fixedsize = true
       width = 1.3]                      
  
A   [label = 'EF1'] 
B1   [label = 'gene by gene\nmeanCt']
B2  [label = 'discard\nall genes',
fontcolor = red] 
C1   [label = 'NA']
C2   [label = '>Detect.Thresh.']
C3   [label = '<Detect.Thresh.']
D1   [label = '2 NA']
D2  [label = '1 NA']
#D3   [label = 'impute',
#fontcolor = red]
D4   [label = 'diff.Ct']
#E1   [label = 'impute',
#fontcolor = red]
E2   [label = 'remaining Ct\n<32 Ct']
E3   [label = 'remaining Ct\n>32 Ct']
E4   [label = 'discard higher\nCt value',
fontcolor = red]
E5   [label = 'valid:\nuse mean Ct',
fontcolor = red]
E6   [label = 'Invalid',
fontcolor = red]
F1   [label = 'discard NA Ct',
fontcolor = red]
#F2   [label = 'impute',
#fontcolor = red]
F3   [label = 'rerun RTqPCR\ndiscard\nuse mean',
fontcolor = darkgreen
fontsize=9]
IMP   [label = 'impute\nrelative conc.',
fontcolor = red]
  
  # edges
  #######
A   -> B1 [label = 'valid:\nboth Ct<mean\u00b12*sd, diff.Ct<0.5',
            style= dashed]
A   -> B2 [label = 'invalid',
                        fontcolor = red]
B1 -> C1
B1 -> C2
B1 -> C3
C1 -> D1
C1 -> D2
#C2 -> D3
C3 -> D4
#D1 -> E1
D2 -> E2
D2 -> E3
E2 -> F1
#E3 -> F2
D4 -> E4 [label = '>3,\nTechnical error']
D4 -> E5 [label = '<Poisson Threshold']
D4 -> E6 [label = 'Poisson threshold < diff.Ct < 3']
E6 -> F3
  
  # grouped edge
  {D1 C2 E3} -> IMP [label = '']
}
")


#############################
