import argparse
import matplotlib.pyplot as plt
from trochoid import trochoid

def draw(
    figsize: list[int],
    output: str,
    shape: str,
    rc: float,
    rm: float,
    rd: float,
    num: int
):
    fig, ax = plt.subplots(
        figsize=(figsize[0]/72, figsize[1]/72)
    )

    outer = (shape == 'hypocycloid')

    x, y = trochoid.ctrochoid(rc, rm, rd, num=num, outer=outer)
    ax.plot(x, y)
    ax.set_aspect('equal', 'box')
    ax.set_title(
        f"{shape}-{rc:0.2f}-{rm:0.2f}-{rd:0.2f}"
    )
    plt.savefig(output)


def main():
    parser = argparse.ArgumentParser()

    fig = parser.add_argument_group('figure')
    fig.add_argument('--figsize', nargs=2, help='figure size w x h',
                     type=int, default=[800, 800])
    fig.add_argument('-o', '--output', default='trochoid.png')

    shape = parser.add_argument_group('shape')
    shape.add_argument('--shape', default='hypocycloid',
                       help='hypocycloid or epicycloid')

    shape.add_argument('--rc', type=float,
                       help='constant circle radius', default=1.0)
    shape.add_argument('--rm', type=float,
                       help='moving circle radius', default=1.0)
    shape.add_argument('--rd', type=float,
                       help='drawing arm length', default=1.0)
    shape.add_argument('-n', '--num', type=int,
                       help='num of points', default=1024)

    args = parser.parse_args()
    draw(
        args.figsize,
        args.output,
        args.shape,
        args.rc,
        args.rm,
        args.rd,
        args.num
    )


if __name__ == '__main__':
    main()
