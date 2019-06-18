bl_info = {
    "name": "Convex Hull Generation in Julia",
    "author": "Jackson Meade (NCSSM) and Jack Snoeyink (UNC Chapel Hill)",
    "version": (1, 0),
    "blender": (2, 80, 0),
    "location": "View3D > Add > Mesh > CGinJulia Convex Hull",
    "description": "Adds a new Mesh Object",
    "warning": "Please make sure you are familiar with the CGinJulia package before using this EXPERIMENTAL add-on",
    "wiki_url": "github.com/snoeyink/CGinJulia",
    "category": "Add Convex Hull",
}


import bpy
from bpy.types import Operator
from bpy.props import FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from mathutils import Vector


def add_object(self, context):
    scale_x = self.scale.x
    scale_y = self.scale.y
    scale_z = self.scale.z

    verts = [Vector((7.137715429456617 * scale_x,46.71163655370694 * scale_y,43.11444831468791 * scale_z)),Vector((59.77347853843603 * scale_x,29.693901997304273 * scale_y,1.4916249683577743 * scale_z)),Vector((45.55810076953275 * scale_x,56.263394323992706 * scale_y,48.06521206481839 * scale_z)),Vector((7.137715429456617 * scale_x,35.78025274030838 * scale_y,46.71163655370694 * scale_z)),Vector((59.77347853843603 * scale_x,57.58037451341368 * scale_y,29.693901997304273 * scale_z)),Vector((45.55810076953275 * scale_x,31.347593001287294 * scale_y,56.263394323992706 * scale_z)),Vector((7.137715429456617 * scale_x,59.43886742717972 * scale_y,35.78025274030838 * scale_z)),Vector((59.77347853843603 * scale_x,56.71235570170625 * scale_y,57.58037451341368 * scale_z)),Vector((45.55810076953275 * scale_x,17.27465626019796 * scale_y,31.347593001287294 * scale_z)),Vector((7.137715429456617 * scale_x,19.53762894376259 * scale_y,10.233634126310234 * scale_z)),Vector((59.77347853843603 * scale_x,30.896467070539344 * scale_y,25.3174853879948 * scale_z)),Vector((45.55810076953275 * scale_x,36.85935066284393 * scale_y,9.899746015333957 * scale_z)),Vector((7.137715429456617 * scale_x,10.233634126310234 * scale_y,6.882215615612899 * scale_z)),Vector((59.77347853843603 * scale_x,25.3174853879948 * scale_y,58.71457895298704 * scale_z)),Vector((45.55810076953275 * scale_x,9.899746015333957 * scale_y,5.29219641246756 * scale_z)),Vector((7.137715429456617 * scale_x,43.11444831468791 * scale_y,19.53762894376259 * scale_z)),Vector((59.77347853843603 * scale_x,1.4916249683577743 * scale_y,30.896467070539344 * scale_z)),Vector((45.55810076953275 * scale_x,48.06521206481839 * scale_y,36.85935066284393 * scale_z)),Vector((46.71163655370694 * scale_x,57.589198037295716 * scale_y,56.396327871218695 * scale_z)),Vector((29.693901997304273 * scale_x,47.07087227611106 * scale_y,15.00685939281916 * scale_z)),Vector((56.263394323992706 * scale_x,21.164049813393092 * scale_y,13.201274395778611 * scale_z)),Vector((46.71163655370694 * scale_x,56.396327871218695 * scale_y,43.11444831468791 * scale_z)),Vector((29.693901997304273 * scale_x,15.00685939281916 * scale_y,1.4916249683577743 * scale_z)),Vector((56.263394323992706 * scale_x,13.201274395778611 * scale_y,48.06521206481839 * scale_z)),Vector((46.71163655370694 * scale_x,43.11444831468791 * scale_y,7.137715429456617 * scale_z)),Vector((29.693901997304273 * scale_x,1.4916249683577743 * scale_y,59.77347853843603 * scale_z)),Vector((56.263394323992706 * scale_x,48.06521206481839 * scale_y,45.55810076953275 * scale_z)),Vector((42.06917206474413 * scale_x,10.233634126310234 * scale_y,36.3513423090245 * scale_z)),Vector((12.371535495951612 * scale_x,25.3174853879948 * scale_y,5.078833173126118 * scale_z)),Vector((12.55018152573494 * scale_x,9.899746015333957 * scale_y,19.73078362505968 * scale_z)),Vector((6.882215615612899 * scale_x,7.137715429456617 * scale_y,10.233634126310234 * scale_z)),Vector((58.71457895298704 * scale_x,59.77347853843603 * scale_y,25.3174853879948 * scale_z)),Vector((5.29219641246756 * scale_x,45.55810076953275 * scale_y,9.899746015333957 * scale_z)),Vector((6.882215615612899 * scale_x,59.43886742717972 * scale_y,7.137715429456617 * scale_z)),Vector((58.71457895298704 * scale_x,56.71235570170625 * scale_y,59.77347853843603 * scale_z)),Vector((5.29219641246756 * scale_x,17.27465626019796 * scale_y,45.55810076953275 * scale_z)),Vector((6.882215615612899 * scale_x,47.832413402500976 * scale_y,59.43886742717972 * scale_z)),Vector((58.71457895298704 * scale_x,26.266415973802854 * scale_y,56.71235570170625 * scale_z)),Vector((5.29219641246756 * scale_x,0.8708468456366525 * scale_y,17.27465626019796 * scale_z)),Vector((6.882215615612899 * scale_x,10.233634126310234 * scale_y,47.832413402500976 * scale_z)),Vector((58.71457895298704 * scale_x,25.3174853879948 * scale_y,26.266415973802854 * scale_z)),Vector((5.29219641246756 * scale_x,9.899746015333957 * scale_y,0.8708468456366525 * scale_z)),Vector((59.43886742717972 * scale_x,47.832413402500976 * scale_y,56.396327871218695 * scale_z)),Vector((56.71235570170625 * scale_x,26.266415973802854 * scale_y,15.00685939281916 * scale_z)),Vector((17.27465626019796 * scale_x,0.8708468456366525 * scale_y,13.201274395778611 * scale_z)),Vector((59.43886742717972 * scale_x,57.589198037295716 * scale_y,46.71163655370694 * scale_z)),Vector((56.71235570170625 * scale_x,47.07087227611106 * scale_y,29.693901997304273 * scale_z)),Vector((17.27465626019796 * scale_x,21.164049813393092 * scale_y,56.263394323992706 * scale_z)),Vector((36.3513423090245 * scale_x,10.233634126310234 * scale_y,43.11444831468791 * scale_z)),Vector((5.078833173126118 * scale_x,25.3174853879948 * scale_y,1.4916249683577743 * scale_z)),Vector((19.73078362505968 * scale_x,9.899746015333957 * scale_y,48.06521206481839 * scale_z)),Vector((36.3513423090245 * scale_x,56.396327871218695 * scale_y,42.06917206474413 * scale_z)),Vector((5.078833173126118 * scale_x,15.00685939281916 * scale_y,12.371535495951612 * scale_z)),Vector((19.73078362505968 * scale_x,13.201274395778611 * scale_y,12.55018152573494 * scale_z)),Vector((47.832413402500976 * scale_x,6.882215615612899 * scale_y,10.233634126310234 * scale_z)),Vector((26.266415973802854 * scale_x,58.71457895298704 * scale_y,25.3174853879948 * scale_z)),Vector((0.8708468456366525 * scale_x,5.29219641246756 * scale_y,9.899746015333957 * scale_z)),Vector((47.832413402500976 * scale_x,59.43886742717972 * scale_y,6.882215615612899 * scale_z)),Vector((26.266415973802854 * scale_x,56.71235570170625 * scale_y,58.71457895298704 * scale_z)),Vector((0.8708468456366525 * scale_x,17.27465626019796 * scale_y,5.29219641246756 * scale_z)),Vector((57.589198037295716 * scale_x,46.71163655370694 * scale_y,59.43886742717972 * scale_z)),Vector((47.07087227611106 * scale_x,29.693901997304273 * scale_y,56.71235570170625 * scale_z)),Vector((21.164049813393092 * scale_x,56.263394323992706 * scale_y,17.27465626019796 * scale_z)),Vector((57.589198037295716 * scale_x,59.43886742717972 * scale_y,56.396327871218695 * scale_z)),Vector((47.07087227611106 * scale_x,56.71235570170625 * scale_y,15.00685939281916 * scale_z)),Vector((21.164049813393092 * scale_x,17.27465626019796 * scale_y,13.201274395778611 * scale_z)),Vector((57.589198037295716 * scale_x,56.396327871218695 * scale_y,46.71163655370694 * scale_z)),Vector((47.07087227611106 * scale_x,15.00685939281916 * scale_y,29.693901997304273 * scale_z)),Vector((21.164049813393092 * scale_x,13.201274395778611 * scale_y,56.263394323992706 * scale_z)),Vector((10.233634126310234 * scale_x,7.137715429456617 * scale_y,19.53762894376259 * scale_z)),Vector((25.3174853879948 * scale_x,59.77347853843603 * scale_y,30.896467070539344 * scale_z)),Vector((9.899746015333957 * scale_x,45.55810076953275 * scale_y,36.85935066284393 * scale_z)),Vector((10.233634126310234 * scale_x,19.53762894376259 * scale_y,43.11444831468791 * scale_z)),Vector((25.3174853879948 * scale_x,30.896467070539344 * scale_y,1.4916249683577743 * scale_z)),Vector((9.899746015333957 * scale_x,36.85935066284393 * scale_y,48.06521206481839 * scale_z)),Vector((10.233634126310234 * scale_x,43.11444831468791 * scale_y,36.3513423090245 * scale_z)),Vector((25.3174853879948 * scale_x,1.4916249683577743 * scale_y,5.078833173126118 * scale_z)),Vector((9.899746015333957 * scale_x,48.06521206481839 * scale_y,19.73078362505968 * scale_z)),Vector((56.396327871218695 * scale_x,46.71163655370694 * scale_y,57.589198037295716 * scale_z)),Vector((15.00685939281916 * scale_x,29.693901997304273 * scale_y,47.07087227611106 * scale_z)),Vector((13.201274395778611 * scale_x,56.263394323992706 * scale_y,21.164049813393092 * scale_z)),Vector((56.396327871218695 * scale_x,42.06917206474413 * scale_y,36.3513423090245 * scale_z)),Vector((15.00685939281916 * scale_x,12.371535495951612 * scale_y,5.078833173126118 * scale_z)),Vector((13.201274395778611 * scale_x,12.55018152573494 * scale_y,19.73078362505968 * scale_z)),Vector((56.396327871218695 * scale_x,59.43886742717972 * scale_y,47.832413402500976 * scale_z)),Vector((15.00685939281916 * scale_x,56.71235570170625 * scale_y,26.266415973802854 * scale_z)),Vector((13.201274395778611 * scale_x,17.27465626019796 * scale_y,0.8708468456366525 * scale_z)),Vector((56.396327871218695 * scale_x,36.3513423090245 * scale_y,43.11444831468791 * scale_z)),Vector((15.00685939281916 * scale_x,5.078833173126118 * scale_y,1.4916249683577743 * scale_z)),Vector((13.201274395778611 * scale_x,19.73078362505968 * scale_y,48.06521206481839 * scale_z)),Vector((56.396327871218695 * scale_x,47.832413402500976 * scale_y,42.06917206474413 * scale_z)),Vector((15.00685939281916 * scale_x,26.266415973802854 * scale_y,12.371535495951612 * scale_z)),Vector((13.201274395778611 * scale_x,0.8708468456366525 * scale_y,12.55018152573494 * scale_z)),Vector((56.396327871218695 * scale_x,43.11444831468791 * scale_y,46.71163655370694 * scale_z)),Vector((15.00685939281916 * scale_x,1.4916249683577743 * scale_y,29.693901997304273 * scale_z)),Vector((13.201274395778611 * scale_x,48.06521206481839 * scale_y,56.263394323992706 * scale_z)),Vector((43.11444831468791 * scale_x,7.137715429456617 * scale_y,46.71163655370694 * scale_z)),Vector((1.4916249683577743 * scale_x,59.77347853843603 * scale_y,29.693901997304273 * scale_z)),Vector((48.06521206481839 * scale_x,45.55810076953275 * scale_y,56.263394323992706 * scale_z)),Vector((43.11444831468791 * scale_x,46.71163655370694 * scale_y,56.396327871218695 * scale_z)),Vector((1.4916249683577743 * scale_x,29.693901997304273 * scale_y,15.00685939281916 * scale_z)),Vector((48.06521206481839 * scale_x,56.263394323992706 * scale_y,13.201274395778611 * scale_z)),Vector((43.11444831468791 * scale_x,56.396327871218695 * scale_y,36.3513423090245 * scale_z)),Vector((1.4916249683577743 * scale_x,15.00685939281916 * scale_y,5.078833173126118 * scale_z)),Vector((48.06521206481839 * scale_x,13.201274395778611 * scale_y,19.73078362505968 * scale_z)),]
    edges = []
    faces = []
    mesh = bpy.data.meshes.new(name="Convex Hull of HelloWorld.py")
    mesh.from_pydata(verts, edges, faces)
    # useful for development when the mesh may be invalid.
    # mesh.validate(verbose=True)
    object_data_add(context, mesh, operator=self)

class OBJECT_OT_add_object(Operator, AddObjectHelper):
    """Generate the Mesh Object from Points and CGinJulia-Constructed Faces"""
    bl_idname = "mesh.add_object"
    bl_label = "Show Representation"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):

        add_object(self, context)

        return {'FINISHED'}


# Registration

def add_object_button(self, context):
    self.layout.operator(
        OBJECT_OT_add_object.bl_idname,
        text="CGinJulia Convex Hull (HelloWorld.py)",
        icon='SURFACE_NSPHERE')

# This allows you to right click on a button and link to the manual
def add_object_manual_map():
    url_manual_prefix = "https://docs.blender.org/manual/en/dev/"
    url_manual_mapping = (
        ("bpy.ops.mesh.add_object", "editors/3dview/object"),
    )
    return url_manual_prefix, url_manual_mapping


def register():
    bpy.utils.register_class(OBJECT_OT_add_object)
    bpy.utils.register_manual_map(add_object_manual_map)
    bpy.types.VIEW3D_MT_mesh_add.append(add_object_button)


def unregister():
    bpy.utils.unregister_class(OBJECT_OT_add_object)
    bpy.utils.unregister_manual_map(add_object_manual_map)
    bpy.types.VIEW3D_MT_mesh_add.remove(add_object_button)


if __name__ == "__main__":
    register()
