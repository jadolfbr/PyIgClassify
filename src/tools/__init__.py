
#Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#This is here because I have written
## #include <utility/vector1.hh>
## utility::vector1< Size >
## A million times.


class vector1(list):
    """
    A list indexed at 1!
    """
    def __init__(self, seq=()):
        list.__init__(self, seq)
        self.insert(0, 0)

    def __len__(self):
        return list.__len__(self) - 1

    def __iter__(self):
        from_one = [x for x in list.__iter__(self)][1:]
        return (x for x in from_one)