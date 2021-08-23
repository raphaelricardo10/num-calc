def func(*args, **kwargs):
    for arg in args:
        print(arg)
    for kw in kwargs:
        print(kw, ":", kwargs[kw])

func("hello", "hi sir", customer="i wanna a coffee, please", sallesman="Good!, Just a second, please")