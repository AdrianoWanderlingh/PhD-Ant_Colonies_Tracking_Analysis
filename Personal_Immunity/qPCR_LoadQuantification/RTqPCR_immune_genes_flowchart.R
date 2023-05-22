pacman::p_load(
  DiagrammeR, # for flow diagrams
  networkD3, # For alluvial/Sankey diagrams
  tidyverse
) # data management and visualization




DiagrammeR::grViz("               # All instructions are within a large character string
digraph immune_genes_data_cleaning_pipeline {    # 'digraph' means 'directional graph', then the graph name 
  
  # graph statement
  #################
  graph [layout = dot,
         splines=line,
         rankdir = TB,            # layout top-to-bottom
         fontsize = 10,
         fontname = 'Times New Roman']
  

  # nodes
  #################
  node [shape = box,           # shape = circle
       fixedsize = true,
       width = 1.3,
      style='filled',
      fillcolor='WhiteSmoke',
      bgcolor='WhiteSmoke',
      color='DimGray'
]                      
  
A   [label = 'Housekeeping\ngene', fillcolor ='Gray70'] 
B1   [label = 'gene by gene\nmean Ct', fillcolor ='Gray70']
B2  [label = 'discard\nall genes',
fontcolor = red] 

C1   [label = 'NA']
C2   [label = ' ≥LOD']
C3   [label = '<LOD']

D1   [label = '2 NA']
D2  [label = '1 NA']
#D3   [label = 'impute',
#fontcolor = red]
D4   [label = 'Diff.Ct', fillcolor ='Gray70']
#E1   [label = 'impute',
#fontcolor = red]
E2   [label = 'remaining Ct\n<(LOD - T.E.)',
#fontsize=9
]
E3   [label = 'remaining Ct\n ≥(LOD - T.E.)',
#fontsize=9
]
E7   [label = '≤Poisson\nThreshold']
E8   [label = '>Poisson\nThreshold']
E5   [label = 'valid:\nuse mean Ct',
fontcolor = red]
E6   [label = 'Invalid\n(discard)',
fontcolor = red]
F1   [label = 'discard NA Ct',
fontcolor = red]
#F2   [label = 'impute',
#fontcolor = red]
IMP   [label = 'impute\nrelative conc.',
fontcolor = red]
J3 [label = '',fillcolor=white,color=white] #empty node used for spacing

#Notes [label = 'Ct: Cycle Threshold; LOD: Limit of Detection; T.E.: Technical Error ', color = white]
  
  # edges
  #######
edge [color = 'DimGray', arrowhead = 'vee']
A   -> B1 [label = 'valid']
A   -> B2 [label = 'invalid']
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
D4 -> E7
E7 -> E5
D4 -> E8
E8 -> E6 
#E6 -> F3
IMP -> J3 [color = white]

  # grouped edge
  {D1 C2 E3} -> IMP [label = '']

#ranking of parts of graph
{rank = same; C1; C2; C3;}
{rank = same; D1; D2; D4;}
{rank = same; E2; E3;}
{rank = same; IMP; J3; F1; E5; E6}
#{rank = sink; Notes; }
}
")

# # valid:\n mean-2sd<Ct<mean+2sd\nΔCt<0.5\nneither NA


## #F3   [label = 'rerun RTqPCR\ndiscard\nuse mean', fontcolor = darkgreen, fontsize=9]

## #E4   [label = 'discard higher\nCt value', fontcolor = red]


#E3 -> F2
#D4 -> E4 [label = '>T.E.']
#D4 -> E6 [label = 'Poisson threshold < ΔCt ≤ T.E.']


#############################
