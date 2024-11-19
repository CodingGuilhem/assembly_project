import click

@click.command()
@click.option('-k','--kmer_size', default=31, help='Kmer size used for the alignment')
def main(count, name):
    """Simple program that greets NAME for a total of COUNT times."""
    for _ in range(count):
        click.echo(f'Hello, {name}!')

if __name__ == '__main__':
    main()