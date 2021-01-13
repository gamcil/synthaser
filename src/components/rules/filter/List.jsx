import { FilterItem } from './Item'

export const FilterList = props => (
  <div>
    <div>
      <button type="button" onClick={props.handleAdd}>Add</button>
    </div>
    <ul>
      {props.filters.map((filter, index) => (
        <FilterItem
          key={filter.uuid}
          data={filter}
          rule={props.rule}
          domains={props.domains}
          handleChange={props.handleChange(index)}
          handleRemove={props.handleRemove(index)}
        />
      ))}
    </ul>
  </div>
)
