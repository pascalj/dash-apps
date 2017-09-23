var param = {
"name": "patviz_params",
"type": "group",
"content": [
{"name":"pattern","label":"Pattern","type":"selection","values":[{"name":"summa","label":"SummaPattern"},{"name":"block","label":"BlockPattern"},{"name":"shift","label":"ShiftTilePattern"},{"name":"seq","label":"SeqTilePattern"}]},
{"name":"numDim","label":"Number of Dimensions","type":"range","value":2,"min":2,"max":2},
 {"type":"br"},
{"name":"blocked_display","label":"Blocked display","type":"checkbox"},
{"name":"balance_extents","label":"Balance extents","type":"checkbox"},
 {"type":"br"},
{"name":"dim_group","type":"group","quantity":"numDim","content":[
 {"name":"size","label":"Extent Dimension","type":"logrange"},
 {"name":"units","label":"Units Dimension","type":"number"},
 {"name":"tile","label":"Tile Dimension","type":"number"}
]}
]
};
