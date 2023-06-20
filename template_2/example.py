class Example:
    def __init__(self, message):
        self.__message = message

    def demo(self):
        print("Reversed message: ", self.__message[::-1])