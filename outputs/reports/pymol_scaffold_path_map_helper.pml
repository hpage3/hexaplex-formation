# PyMOL scaffold path map helper
# Load this script in PyMOL, then compare these generated path colors
# against the known colored strand/path representation.
# This helper validates visual correspondence only; it does not prove biological truth.
load outputs/intermediates/normalized_structures/hexaplex_scaffold_only_complement_heavy_deduped.pdb, scaffold_path_map_model
hide everything, scaffold_path_map_model
show cartoon, scaffold_path_map_model
show sticks, scaffold_path_map_model
color gray70, scaffold_path_map_model
# Map CSV: inputs/metadata/scaffold_path_map_candidate.csv

# strand_id=1 strand_label=block_1 residues=30
select path_block_1, scaffold_path_map_model and ((resi 1+2+3+4+5+6+7+8+9+10+11+12+13+14+15+16+17+18+19+20+21+22+23+24+25+26+27+28+29+30))
color red, path_block_1
show sticks, path_block_1
label first path_block_1 and name CA, "block_1"
print "block_1: 30 residue(s), selection path_block_1"

# strand_id=2 strand_label=block_2 residues=30
select path_block_2, scaffold_path_map_model and ((resi 31+32+33+34+35+36+37+38+39+40+41+42+43+44+45+46+47+48+49+50+51+52+53+54+55+56+57+58+59+60))
color orange, path_block_2
show sticks, path_block_2
label first path_block_2 and name CA, "block_2"
print "block_2: 30 residue(s), selection path_block_2"

# strand_id=3 strand_label=block_3 residues=30
select path_block_3, scaffold_path_map_model and ((resi 61+62+63+64+65+66+67+68+69+70+71+72+73+74+75+76+77+78+79+80+81+82+83+84+85+86+87+88+89+90))
color yellow, path_block_3
show sticks, path_block_3
label first path_block_3 and name CA, "block_3"
print "block_3: 30 residue(s), selection path_block_3"

# strand_id=4 strand_label=block_4 residues=30
select path_block_4, scaffold_path_map_model and ((resi 91+92+93+94+95+96+97+98+99+100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118+119+120))
color green, path_block_4
show sticks, path_block_4
label first path_block_4 and name CA, "block_4"
print "block_4: 30 residue(s), selection path_block_4"

# strand_id=5 strand_label=block_5 residues=30
select path_block_5, scaffold_path_map_model and ((resi 121+122+123+124+125+126+127+128+129+130+131+132+133+134+135+136+137+138+139+140+141+142+143+144+145+146+147+148+149+150))
color cyan, path_block_5
show sticks, path_block_5
label first path_block_5 and name CA, "block_5"
print "block_5: 30 residue(s), selection path_block_5"

# strand_id=6 strand_label=block_6 residues=30
select path_block_6, scaffold_path_map_model and ((resi 151+152+153+154+155+156+157+158+159+160+161+162+163+164+165+166+167+168+169+170+171+172+173+174+175+176+177+178+179+180))
color blue, path_block_6
show sticks, path_block_6
label first path_block_6 and name CA, "block_6"
print "block_6: 30 residue(s), selection path_block_6"

orient scaffold_path_map_model
zoom scaffold_path_map_model
# If generated colors do not match the reference colored paths, edit
# inputs/metadata/scaffold_path_map_manual_template.csv into
# inputs/metadata/scaffold_path_map_manual.csv and rerun the workflow.
