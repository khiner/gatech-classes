# Notes for 8803

## Inspiration

A collection of inspiring GNN-related papers:

(Also see [mesh2audio Notes.md](https://github.com/khiner/mesh2audio/blob/main/Notes.md) - this ties in since it's in mesh land!)

- Learning Mesh-Based Simulation with Graph Networks
  - Paper: https://arxiv.org/abs/2010.03409
  - Video: https://youtu.be/KfZFgSff9N8?si=6DxM6LYEDGSRh2RI
  - This one _is_ GNN-based.
- Learning to Simulate Complex Physics with Graph Networks
  - In the MLG reading list.
- Temporal Graph Networks for Deep Learning on Dynamic Graphs
  - In the MLG reading list.

## GeoLDM

I have the "hacker" role for this paper.

I ran this on my MBP and it took about a day to do (on the CPU):

```shell
$ python eval_sample.py --model_path outputs/qm9_latent2 --n_samples 100
```

Outputs are in `GeoLDM/GeoLDM/outputs/qm9_latent2/eval_sep_9_2023/`.
See `GeoLDMRunOutput.txt` for the full output.

So at the very least, I have locally run gifs to share, and I can make a 1-cell jupyter notebook that just runs the original model.

Of course, I'll need to do more than this.
Not sure what yet.

The dream is real-time learning/visualization.
So essentially, imagine any of the output gifs from `eval_sample.py` learned on-the-fly.

**Goal: Run GeoLDM on Metal.**

My ~1-day-long run for 100 samples was CPU-only.
Let's see how much we can speed it up still on Mac!

- Get PyTorch Metal acceleration working:
  - https://pytorch.org/docs/stable/notes/mps.html
  - Installed nightly with
  `pip3 install --pre torch torchvision torchaudio --index-url https://download.pytorch.org/whl/nightly/cpu`

Timing with mps:

Sampling handful of molecules.
Time:  164.87406301498413
Sampling stable molecules.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Time:  285.89433693885803

Sampling handful of molecules.
Time:  162.48754715919495
Sampling stable molecules.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Time:  285.4007008075714

Timing with cpu:

Sampling handful of molecules.
Time:  169.80601906776428
Sampling stable molecules.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Found stable mol.
Time:  336.21314120292664
