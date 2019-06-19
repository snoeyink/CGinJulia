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

    verts = [Vector((1.2237668 * scale_x,2.2781324 * scale_y,17.893223 * scale_z)),Vector((45.35256 * scale_x,33.612114 * scale_y,50.394894 * scale_z)),Vector((40.321907 * scale_x,7.8105497 * scale_y,1.0372853 * scale_z)),Vector((58.25678 * scale_x,7.150562 * scale_y,59.202316 * scale_z)),Vector((13.813097 * scale_x,0.13966084 * scale_y,3.514216 * scale_z)),Vector((4.352095 * scale_x,36.28717 * scale_y,3.550601 * scale_z)),Vector((31.95387 * scale_x,23.363207 * scale_y,3.7901187 * scale_z)),Vector((0.9114218 * scale_x,0.53810835 * scale_y,23.800407 * scale_z)),Vector((41.59571 * scale_x,52.812653 * scale_y,56.535988 * scale_z)),Vector((16.72496 * scale_x,4.137032 * scale_y,26.386692 * scale_z)),Vector((23.11836 * scale_x,53.925125 * scale_y,53.223324 * scale_z)),Vector((28.068716 * scale_x,43.18795 * scale_y,7.7776837 * scale_z)),Vector((54.739437 * scale_x,2.58497 * scale_y,31.99495 * scale_z)),Vector((50.802284 * scale_x,53.71783 * scale_y,5.2884007 * scale_z)),Vector((6.441815 * scale_x,3.6299586 * scale_y,51.589935 * scale_z)),Vector((47.477272 * scale_x,26.334743 * scale_y,25.156775 * scale_z)),Vector((44.2535 * scale_x,20.764626 * scale_y,20.085047 * scale_z)),Vector((41.31904 * scale_x,28.373695 * scale_y,25.169785 * scale_z)),Vector((47.932526 * scale_x,23.969786 * scale_y,29.911337 * scale_z)),Vector((46.32869 * scale_x,13.890867 * scale_y,1.7579842 * scale_z)),Vector((49.120247 * scale_x,13.216646 * scale_y,35.849182 * scale_z)),Vector((37.203117 * scale_x,59.44161 * scale_y,48.391487 * scale_z)),Vector((8.837421 * scale_x,23.998074 * scale_y,11.243771 * scale_z)),]
    edges = []
    faces = [[1,6,8],[1,8,5],[4,13,15],[4,14,13],[5,1,8],[5,6,1],[5,13,3],[6,1,5],[6,5,3],[6,8,1],[6,11,15],[6,15,8],[9,4,15],[9,14,4],[13,3,5],[13,14,20],[13,15,4],[13,20,3],[14,6,3],[14,9,22],[14,20,13],[15,6,11],[15,8,6],[20,3,13],[22,9,11],[22,11,6],[22,14,9],]
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
