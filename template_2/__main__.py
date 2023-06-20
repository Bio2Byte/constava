import argparse
from template_2.example import Example


class TemplateMain:
    def __init__(self):
        pass

    def run(self, parsed_args):
        Example(parsed_args.message).demo()

def main():
    print("Constava 0.0.1 - Template 2")
    parser = argparse.ArgumentParser(description='Constava - Example')
    parser.add_argument('-m', '--message', default='hello world')

    parsed_args = parser.parse_args()
    TemplateMain().run(parsed_args)

if __name__ == "__main__":
    main()
