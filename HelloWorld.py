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

    verts = [Vector((84 * scale_x,43 * scale_y,70 * scale_z)),Vector((23 * scale_x,42 * scale_y,30 * scale_z)),Vector((57 * scale_x,45 * scale_y,45 * scale_z)),Vector((44 * scale_x,64 * scale_y,76 * scale_z)),Vector((90 * scale_x,67 * scale_y,64 * scale_z)),Vector((63 * scale_x,88 * scale_y,23 * scale_z)),Vector((64 * scale_x,87 * scale_y,54 * scale_z)),Vector((7 * scale_x,75 * scale_y,63 * scale_z)),Vector((13 * scale_x,67 * scale_y,80 * scale_z)),Vector((97 * scale_x,7 * scale_y,50 * scale_z)),Vector((78 * scale_x,43 * scale_y,12 * scale_z)),Vector((38 * scale_x,62 * scale_y,8 * scale_z)),Vector((82 * scale_x,76 * scale_y,13 * scale_z)),Vector((34 * scale_x,27 * scale_y,15 * scale_z)),Vector((22 * scale_x,33 * scale_y,67 * scale_z)),Vector((99 * scale_x,41 * scale_y,98 * scale_z)),Vector((51 * scale_x,19 * scale_y,63 * scale_z)),Vector((88 * scale_x,10 * scale_y,42 * scale_z)),Vector((55 * scale_x,49 * scale_y,60 * scale_z)),Vector((92 * scale_x,57 * scale_y,12 * scale_z)),Vector((88 * scale_x,50 * scale_y,34 * scale_z)),Vector((30 * scale_x,23 * scale_y,54 * scale_z)),Vector((100 * scale_x,24 * scale_y,73 * scale_z)),]
    edges = []
    faces = [[1,7,14],[1,11,7],[1,13,11],[1,14,13],[4,6,12],[4,12,19],[4,15,6],[4,19,15],[5,6,7],[5,7,11],[5,11,12],[5,12,6],[6,8,7],[6,15,8],[7,8,14],[8,15,14],[9,10,17],[9,15,22],[9,16,15],[9,17,21],[9,19,10],[9,21,16],[9,22,19],[10,11,13],[10,13,17],[10,19,11],[11,19,12],[13,14,21],[13,21,17],[14,15,16],[14,16,21],[15,19,22],]
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
