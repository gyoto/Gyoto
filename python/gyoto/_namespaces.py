'''Namespace wrapper helpers

Not part of the API

'''

def make_namespace(Generic, _globals):
    import sys as _sys
    import inspect as _inspect
    from gyoto.core import requirePlugin as _requirePlugin, havePlugin as _havePlugin
    import gyoto.core as _core
    import gyoto.std as _std
    
    _modules = [_core, _std]
    
    _requirePlugin("lorene", True)
    if _havePlugin("lorene"):
        import gyoto.lorene as _lorene
        _modules.append(_lorene)
        
    res = []
    
    for _mod in _modules:
        for _name, _obj in _inspect.getmembers(_mod):
            if _inspect.isclass(_obj) and issubclass(_obj, Generic):
                __import__(_obj.__module__, globals(), locals(), [_obj.__name__,])
                _globals.update({_obj.__name__:
                                 getattr(_sys.modules[_obj.__module__], 
                                         _obj.__name__)})
                res . append(_obj.__name__)
    
    
    return res
