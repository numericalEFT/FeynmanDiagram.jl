# Parquet Algorithm to Build Diagrams

## API

```@autodocs
Modules = [FeynmanDiagram.Parquet]
```

## Usage

```@repl
using FeynmanDiagram
para = DiagPara(type = Ver4Diag, innerLoopNum = 1, hasTau = true);
Parquet.vertex4(para)

para = DiagPara(type = Ver3Diag, innerLoopNum = 1, hasTau = true);
Parquet.vertex3(para)

para = DiagPara(type = SigmaDiag, innerLoopNum = 1, hasTau = true);
Parquet.sigma(para)

para = DiagPara(type = PolarDiag, innerLoopNum = 1, hasTau = true);
Parquet.polarization(para)
```