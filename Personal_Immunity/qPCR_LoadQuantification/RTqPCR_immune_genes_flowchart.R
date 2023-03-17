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
  node [shape = circle,           # shape = circle
       fixedsize = true
       width = 1.3]                      
  
A   [label = 'EF1'] 
B1   [label = 'meanCt']
B2  [label = 'discard\nall genes',
fontcolor = red] 
C1   [label = 'NA']
C2   [label = '>35 Ct']
C3   [label = '<35 Ct']
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
fontcolor = green]
IMP   [label = 'impute',
fontcolor = red]
  
  # edges
  #######
A   -> B1 [label = 'valid:\nboth Ct<mean\u00b1sd, diff.Ct<0.5']
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





##############################


DiagrammeR::grViz('
digraph surveillance_diagram {
  
  # graph statement
  #################
  graph [layout = dot,
         rankdir = TB,
         fontsize = 10]
  
  # nodes
  #######
  node [shape = circle,
        fixedsize = true,
        width = 1.3]
  
  A [label = "EF1"]
  B1 [label = "meanCt"]
  B2 [label = "discard\nall genes",
      fontcolor = red]
  C1 [label = "NA"]
  C2 [label = ">35 Ct"]
  C3 [label = "<35 Ct"]
  D1 [label = "2 NA"]
  D2 [label = "1 NA"]
  D4 [label = "diff Ct"]
  E2 [label = "remaining Ct\n<32 Ct"]
  E3 [label = "remaining Ct\n>32 Ct"]
  E4 [label = "discard higher\nCt value",
      fontcolor = red]
  E5 [label = "valid:\nuse mean Ct",
      fontcolor = green]
  E6 [label = "Invalid",
      fontcolor = red]
  F1 [label = "discard NA Ct",
      fontcolor = red]
  F3 [label = "rerun RTqPCR\ndiscard\nuse mean",
      fontcolor = green]
  IMP [label = "impute",
       fontcolor = red]
  
  # edges
  #######
  A -> B1 [label = "valid:\nboth Ct<mean+-sd, diff Ct<0.5"]
  A -> B2 [label = "invalid",
            fontcolor = red]
  B1 -> C1
  B1 -> C2
  B1 -> C3
  C1 -> D1
  C1 -> D2
  D2 -> E2
  D2 -> E3
  E2 -> F1
  D4 -> E4 [label = ">3,\nTechnical error"]
  D4 -> E5 [label = "<Poisson Threshold"]
  D4 -> E6 [label = "Poisson threshold < diff Ct < 3"]
  E6 -> F3
  
  # nodes and edges to be highlighted
  ###################################
  {D1 C2 E3} -> IMP [label = "",
                     fontcolor = darkgreen,
                     color = darkgreen,
                     style = dashed]
  
  # global node and edge styles
  #############################
  node [fontname = "Helvetica",
        fontsize = 9,
        shape = circle,
        style = filled,
        fillcolor = white,
        color = black,
        penwidth = 1.0]
  edge [fontname = "Helvetica",
        fontsize = 8,
        color = black,
        penwidth = 1.0,
        arrowsize = 0.7]
}
')
