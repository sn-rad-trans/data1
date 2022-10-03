#! /bin/csh

### inputs at 2d with default settings (see script)
python mk_snia_toy_model.py --highni --noplot
python mk_snia_toy_model.py --lowni --noplot

### inputs at 1h used in STELLA runs
python mk_snia_toy_model.py --mtot 1.0 --ekin 1.0 --dvel 200.0 --densprof 'expon' --tend 4.1667e-2 \
    --tempmin 5e3 --mige 0.0 --dmfine 1e-4 --transprof 'cosine' --xfracti 0.0 --ime 'ca,s,si' --xfracime '0.1,0.35,0.55' \
    --mni56 0.6 --dmni56 0.4 --mime 0.4 --fout 'snia_toy06_1h_lowres.dat' --noplot

python mk_snia_toy_model.py --mtot 1.0 --ekin 1.0 --dvel 200.0 --densprof 'expon' --tend 4.1667e-2 \
    --tempmin 5e3 --mige 0.0 --dmfine 1e-4 --transprof 'cosine' --xfracti 0.0 --ime 'ca,s,si' --xfracime '0.1,0.35,0.55' \
    --mni56 0.1 --dmni56 0.2 --mime 0.9 --fout 'snia_toy01_1h_lowres.dat' --noplot

exit
