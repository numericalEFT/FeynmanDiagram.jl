# Parquet Algorithm to Build Diagrams

## API

```@autodocs
Modules = [FeynmanDiagram.Parquet]
```

## Usage

```@repl
using FeynmanDiagram
para = DiagPara(diagType = Ver4Diag, innerLoopNum = 1, hasTau = true);
Parquet.vertex4(para)

para = DiagPara(diagType = Ver3Diag, innerLoopNum = 1, hasTau = true);
Parquet.vertex3(para)

para = DiagPara(diagType = SigmaDiag, innerLoopNum = 1, hasTau = true);
Parquet.sigma(para)

para = DiagPara(diagType = PolarDiag, innerLoopNum = 1, hasTau = true);
Parquet.polarization(para)
```