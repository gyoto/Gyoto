# Helpers for the Gyoto test suite
import inspect
from gyoto.core import Factory, Error

def test_xmlio(self, namespace, invalidator):
    for classname, cls in inspect.getmembers(namespace):
        if (invalidator(classname, cls)):
            continue
        with self.subTest(f"{classname=}"):
            try:
                obj=cls()
            except Error as e:
                raise Exception(
                    f'unable to instanciate {classname}: \n' +
                    e.get_message()
                )
            try:
                a = str(obj)
            except Error as e:
                raise Exception(
                    f'unable to describe {classname} as XML: \n' +
                    e.get_message()
                )
            try:
                f = Factory(a)
                factory_getter_name = f.kind().lower()
                obj2 = getattr(f, factory_getter_name)()
            except Error as e:
                raise Exception(
                    f'unable to construct {classname} from XML: \n' +
                    e.get_message()
                )
            try:
                b = str(obj2)
            except Error as e:
                raise Exception(
                    f'unable to describe second instance as XML: \n' +
                    e.get_message()
                )
            self.assertEqual(
                a, b,
                f"XML roudtrip is not idempotent for class {classname}"
            )
            del obj
            del obj2

def test_exposed_classes(self, plg, module, nspace, special):
    items = nspace.Generic.registeredPluginsSlashKinds()
    classes = [item.split('/')[1] for item in items
               if item.startswith(plg + '/')]
    for classname in classes:
        try:
            expname = special[classname]
        except KeyError:
            expname = classname
        self.assertTrue(
            hasattr(module, expname)
            and issubclass(getattr(module, expname), nspace.Generic),
            repr(module) + "does not expose class " + classname
        )
        self.assertTrue(
            hasattr(nspace, classname)
            and issubclass(getattr(nspace, classname), nspace.Generic),
            repr(nspace) + "does not expose class " + classname
        )
        self.assertEqual(
            getattr(module, expname), getattr(nspace, classname),
            f"class {classname} is not consistently exposed in "
            f"{repr(module)} and in {repr(nspace)}"
        )
