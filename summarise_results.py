


# Python implementation to
# read last N lines of a file
# through Exponential search

# Function to read
# last N lines of the file
def LastNlines(fname, N):

    # assert statement check
    # a condition
    assert N >= 0

    # declaring variable
    # to implement
    # exponential search
    pos = N + 1

    # list to store
    # last N lines
    lines = []

    # opening file using with() method
    # so that file get closed
    # after completing work
    with open(fname) as f:

        # loop which runs
        # until size of list
        # becomes equal to N
        while len(lines) <= N:

            # try block
            try:
                # moving cursor from
                # left side to
                # pos line from end
                f.seek(-pos, 2)

            # exception block
            # to handle any run
            # time error
            except IOError:
                f.seek(0)
                break

            # finally block
            # to add lines
            # to list after
            # each iteration
            finally:
                lines = list(f)

            # increasing value
            # of variable
            # exponentially
            pos *= 2

    # returning the
    # whole list
    # which stores last
    # N lines
    return lines[-N:]


import glob
files = glob.glob("./*.out")
data = [LastNlines(filename,10) for filename in files]
#print(data)

def catch(func, handle=lambda e : e, *args, **kwargs):
    try:
        return func(*args, **kwargs)
    except Exception as e:
        return handle(e)

splits = [catch(lambda : x[-1].split("]")[-2].split("'")[-1].split(",")) for x in data]
print(splits)

metadata = [catch(lambda : (x[1],x[2],x[3])) for x in splits]

print(metadata)

d={}
for x in metadata:
    d[x] = d.get(x, 0) + 1

e = {}
for x in d:
    if isinstance(x,tuple):
        e[x] = d[x]

total = sum(e.values())
for x in e:
    e[x] /= total

print(e)

mylist = list(e.items())
mylist.sort(key=lambda a: a[1],reverse=True)
print(mylist)

GD = sum([e[x] for x in e if int(x[1]) >=0 and int(x[2]) == -1])
print(GD)

