from mscg.cli import cgfm

cgfm.main(
    top="cg_topol.top",
    traj="CG.lammpstrj",
    cut=12.0,
    pair=[
        "model=BSpline,type=B1:B1,min=3.0,max=12.0,resolution=0.1",
        "model=BSpline,type=B2:B2,min=3.0,max=12.0,resolution=0.1",
        "model=BSpline,type=B3:B3,min=3.0,max=12.0,resolution=0.1",
        "model=BSpline,type=B1:B2,min=3.0,max=12.0,resolution=0.1",
        "model=BSpline,type=B2:B3,min=3.0,max=12.0,resolution=0.1",
        "model=BSpline,type=B1:SOL,min=3.0,max=12.0,resolution=0.1",
        "model=BSpline,type=B2:SOL,min=3.0,max=12.0,resolution=0.1",
        "model=BSpline,type=B3:SOL,min=3.0,max=12.0,resolution=0.1",
    ],
    bond=[
        "model=BSpline,type=B1:B2,min=1.2,max=1.5,resolution=0.01",
        "model=BSpline,type=B2:B3,min=1.2,max=1.5,resolution=0.01",
    ],
)
