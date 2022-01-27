# Parquet Algorithm to Build Diagrams

## API

```@autodocs
Modules = [FeynmanDiagram.Parquet]
```

## Usage

```@repl
using FeynmanDiagram
para = GenericPara(diagType = Ver4Diag, innerLoopNum = 1, hasTau = true);
Parquet.vertex4(para)

para = GenericPara(diagType = Ver3Diag, innerLoopNum = 1, hasTau = true);
Parquet.vertex3(para)

para = GenericPara(diagType = SigmaDiag, innerLoopNum = 1, hasTau = true);
Parquet.sigma(para)

para = GenericPara(diagType = PolarDiag, innerLoopNum = 1, hasTau = true);
Parquet.polarization(para)
```